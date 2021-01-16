% Calculate joint inertia matrix for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP10_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:14:52
% EndTime: 2021-01-15 19:14:54
% DurationCPUTime: 0.32s
% Computational Cost: add. (452->110), mult. (849->148), div. (0->0), fcn. (869->6), ass. (0->51)
t104 = cos(qJ(3));
t74 = sin(pkin(8));
t75 = cos(pkin(8));
t77 = sin(qJ(3));
t62 = t104 * t74 + t77 * t75;
t107 = 0.2e1 * t62;
t68 = -t75 * pkin(2) - pkin(1);
t106 = 0.2e1 * t68;
t61 = -t104 * t75 + t77 * t74;
t105 = t61 * pkin(4);
t102 = pkin(6) + qJ(2);
t63 = t102 * t74;
t64 = t102 * t75;
t59 = t104 * t64 - t77 * t63;
t78 = cos(qJ(4));
t103 = t78 * t59;
t101 = -qJ(5) - pkin(7);
t100 = t74 ^ 2 + t75 ^ 2;
t76 = sin(qJ(4));
t72 = t76 ^ 2;
t73 = t78 ^ 2;
t99 = t72 + t73;
t98 = MDP(24) * pkin(4);
t97 = qJ(5) * t62;
t96 = t75 * MDP(4);
t95 = t61 * MDP(18);
t66 = t101 * t78;
t94 = t66 * MDP(22);
t69 = -t78 * pkin(4) - pkin(3);
t93 = t69 * MDP(24);
t92 = t76 * MDP(17);
t91 = t76 * MDP(22);
t90 = t78 * MDP(21);
t89 = t76 * t78 * MDP(15);
t57 = t61 * pkin(3) - t62 * pkin(7) + t68;
t53 = t78 * t57 - t76 * t59;
t88 = -MDP(23) * pkin(4) + MDP(16);
t87 = MDP(21) + t98;
t86 = (-MDP(20) - MDP(22)) * t76;
t51 = -t78 * t97 + t105 + t53;
t52 = t103 + (t57 - t97) * t76;
t85 = t51 * t78 + t52 * t76;
t84 = t78 * MDP(19) - t76 * MDP(20);
t83 = MDP(19) * t76 + MDP(20) * t78;
t82 = t76 * MDP(21) + t78 * MDP(22);
t58 = t104 * t63 + t77 * t64;
t81 = t53 * MDP(19) - (t76 * t57 + t103) * MDP(20) - t52 * MDP(22);
t80 = -t90 + t91 + t93;
t65 = t101 * t76;
t55 = t76 * t62 * pkin(4) + t58;
t1 = [MDP(1) + 0.2e1 * pkin(1) * t96 + (t100 * qJ(2) ^ 2 + pkin(1) ^ 2) * MDP(6) + (t51 ^ 2 + t52 ^ 2 + t55 ^ 2) * MDP(24) + (MDP(12) * t106 + t95 + (t78 * MDP(16) - MDP(8) - t92) * t107) * t61 + 0.2e1 * (t51 * MDP(21) + t81) * t61 + (-t85 * MDP(23) + t82 * t55 + t83 * t58) * t107 + 0.2e1 * t100 * MDP(5) * qJ(2) + (MDP(13) * t106 + (t73 * MDP(14) + MDP(7) - 0.2e1 * t89) * t62) * t62; -t96 - pkin(1) * MDP(6) + t85 * MDP(24) + (-t99 * MDP(23) + MDP(13)) * t62 + (MDP(12) + (MDP(19) + MDP(21)) * t78 + t86) * t61; t99 * MDP(24) + MDP(6); -t59 * MDP(13) + (-t51 * t76 + t52 * t78) * MDP(23) + (t51 * t65 - t52 * t66) * MDP(24) + t80 * t55 + (-MDP(12) - t84) * t58 + (t76 * MDP(16) + t78 * MDP(17) + t65 * MDP(21) - t83 * pkin(7) - MDP(10) + t94) * t61 + (MDP(9) + (-t72 + t73) * MDP(15) + (-pkin(3) * MDP(20) + t69 * MDP(22) - t65 * MDP(23)) * t78 + (t78 * MDP(14) - pkin(3) * MDP(19) + t69 * MDP(21) + t66 * MDP(23)) * t76) * t62; (t78 * t65 - t76 * t66) * MDP(24); MDP(11) + t72 * MDP(14) + 0.2e1 * t89 + 0.2e1 * (-t65 * t76 - t66 * t78) * MDP(23) + (t65 ^ 2 + t66 ^ 2) * MDP(24) + (-0.2e1 * t90 + 0.2e1 * t91 + t93) * t69 + 0.2e1 * t84 * pkin(3); t95 + (t53 + 0.2e1 * t105) * MDP(21) + t51 * t98 + (-t92 + (-MDP(21) * qJ(5) + t88) * t78) * t62 + t81; t86 + (MDP(19) + t87) * t78; t94 + (-MDP(20) * pkin(7) + MDP(17)) * t78 + t87 * t65 + (-MDP(19) * pkin(7) + t88) * t76; MDP(18) + (0.2e1 * MDP(21) + t98) * pkin(4); t55 * MDP(24) + t82 * t62; 0; t80; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
