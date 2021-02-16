% Calculate joint inertia matrix for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP12_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:25:55
% EndTime: 2021-01-15 19:25:57
% DurationCPUTime: 0.35s
% Computational Cost: add. (259->112), mult. (455->141), div. (0->0), fcn. (352->4), ass. (0->45)
t63 = cos(qJ(4));
t76 = MDP(20) + MDP(22);
t93 = t76 * t63;
t64 = cos(qJ(3));
t92 = 0.2e1 * t64;
t91 = (pkin(1) * MDP(6));
t61 = sin(qJ(4));
t90 = t61 * t63;
t65 = -pkin(1) - pkin(6);
t89 = t61 * t65;
t88 = t63 * t65;
t87 = -qJ(5) - pkin(7);
t57 = t61 ^ 2;
t59 = t63 ^ 2;
t86 = t57 + t59;
t84 = MDP(24) * pkin(4);
t83 = qJ(5) * t64;
t55 = t87 * t63;
t82 = t55 * MDP(22);
t56 = -t63 * pkin(4) - pkin(3);
t81 = t56 * MDP(24);
t80 = t61 * MDP(17);
t79 = t61 * MDP(22);
t78 = t63 * MDP(21);
t77 = MDP(19) + MDP(21);
t62 = sin(qJ(3));
t75 = t62 * t88;
t74 = MDP(15) * t90;
t73 = -MDP(23) * pkin(4) + MDP(16);
t72 = MDP(21) + t84;
t53 = t62 * pkin(3) - t64 * pkin(7) + qJ(2);
t49 = t63 * t53;
t45 = -t63 * t83 + t49 + (pkin(4) - t89) * t62;
t46 = t75 + (t53 - t83) * t61;
t71 = -t45 * t61 + t46 * t63;
t54 = t87 * t61;
t70 = -t54 * t61 - t55 * t63;
t69 = -MDP(19) * t61 - MDP(20) * t63;
t68 = t61 * MDP(21) + t63 * MDP(22);
t67 = (-t62 * t89 + t49) * MDP(19) - (t61 * t53 + t75) * MDP(20) - t46 * MDP(22);
t66 = -t78 + t79 + t81;
t60 = t64 ^ 2;
t58 = t62 ^ 2;
t50 = (pkin(4) * t61 - t65) * t64;
t1 = [MDP(1) + t58 * MDP(18) + (t45 ^ 2 + t46 ^ 2 + t50 ^ 2) * MDP(24) + ((-2 * MDP(4) + t91) * pkin(1)) + (MDP(13) * t92 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + ((-t45 * t63 - t46 * t61) * MDP(23) + t68 * t50) * t92 + (t59 * MDP(14) + 0.2e1 * t69 * t65 + MDP(7) - 0.2e1 * t74) * t60 + 0.2e1 * (qJ(2) * MDP(12) + (t63 * MDP(16) - MDP(8) - t80) * t64 + t45 * MDP(21) + t67) * t62; MDP(4) - t91 + (-t50 * t64 + t71 * t62) * MDP(24) + (t77 * t61 + t93) * (-t58 - t60); MDP(6) + (t86 * t58 + t60) * MDP(24); t71 * MDP(23) + (t45 * t54 - t46 * t55) * MDP(24) + t66 * t50 + (-t65 * MDP(13) + t61 * MDP(16) + t63 * MDP(17) + t54 * MDP(21) + t69 * pkin(7) - MDP(10) + t82) * t62 + (MDP(9) + t65 * MDP(12) + MDP(14) * t90 + (-t57 + t59) * MDP(15) + (-pkin(3) * t61 + t88) * MDP(19) + (-pkin(3) * t63 - t89) * MDP(20) + (-t54 * t63 + t55 * t61) * MDP(23) + t68 * t56) * t64; (t86 * MDP(23) + t70 * MDP(24) - MDP(13)) * t62 + (-t76 * t61 + t77 * t63 + MDP(12) - t81) * t64; MDP(11) + t57 * MDP(14) + 0.2e1 * t74 + 0.2e1 * t70 * MDP(23) + (t54 ^ 2 + t55 ^ 2) * MDP(24) + (-0.2e1 * t78 + 0.2e1 * t79 + t81) * t56 + 0.2e1 * (t63 * MDP(19) - t61 * MDP(20)) * pkin(3); t45 * t84 + t49 * MDP(21) + (MDP(18) + (0.2e1 * pkin(4) - t89) * MDP(21)) * t62 + (-t80 + (-MDP(21) * qJ(5) + t73) * t63) * t64 + t67; (-t93 + (-MDP(19) - t72) * t61) * t62; t82 + (-MDP(20) * pkin(7) + MDP(17)) * t63 + t72 * t54 + (-MDP(19) * pkin(7) + t73) * t61; MDP(18) + (0.2e1 * MDP(21) + t84) * pkin(4); t50 * MDP(24) + t68 * t64; -t64 * MDP(24); t66; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
