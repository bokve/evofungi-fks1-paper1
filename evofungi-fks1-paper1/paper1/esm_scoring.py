"""Masked marginal primitive; model weights are never loaded by this module."""
from .esm_contract import CANONICAL_AA, window


def masked_tokens(alphabet, sequence, position_1based):
    w = window(sequence, position_1based)
    if not alphabet.prepend_bos or not alphabet.append_eos:
        raise ValueError('ESM-1v/ESM-1b BOS/EOS contract mismatch')
    _, strings, tokens = alphabet.get_batch_converter()([('frozen_FKS1_window', w['sequence'])])
    if strings != [w['sequence']] or tokens.shape != (1, len(w['sequence']) + 2) or tokens.shape[1] > 1024:
        raise ValueError('Tokenizer truncated/reformatted sequence or exceeded token limit')
    idx = w['local_index_0based'] + 1
    wt = sequence[position_1based - 1]
    if wt not in CANONICAL_AA or int(tokens[0,idx]) != alphabet.get_idx(wt) or alphabet.get_idx(wt) == alphabet.unk_idx:
        raise ValueError('WT token/index identity mismatch')
    tokens[0,idx] = alphabet.mask_idx
    if int((tokens == alphabet.mask_idx).sum()) != 1:
        raise ValueError('Exactly one target token must be masked')
    return tokens, idx, wt, w


def score_masked_position(model, alphabet, sequence, position_1based, device='cuda'):
    import torch
    tokens, idx, wt, w = masked_tokens(alphabet, sequence, position_1based)
    if model.training:
        raise ValueError('Scoring requires eval mode')
    if not hasattr(model, 'args') or not hasattr(model.args, 'max_positions') or tokens.shape[1] > model.args.max_positions:
        raise ValueError('Checkpoint token capacity is absent/incompatible; explicit technical review required')
    with torch.inference_mode():
        logits = model(tokens.to(device))['logits']
        if logits.shape[:2] != tokens.shape or logits.shape[-1] != len(alphabet):
            raise ValueError('Checkpoint logits/token shape mismatch')
        logp = torch.log_softmax(logits[0,idx].float(), dim=-1)
        values = {aa: float((logp[alphabet.get_idx(aa)] - logp[alphabet.get_idx(wt)]).cpu()) for aa in sorted(CANONICAL_AA)}
    if not all(torch.isfinite(torch.tensor(v)).item() for v in values.values()) or values[wt] != 0:
        raise ValueError('Nonfinite masked marginal or failed WT self-zero check')
    return values, {key:value for key,value in w.items() if key != 'sequence'}
