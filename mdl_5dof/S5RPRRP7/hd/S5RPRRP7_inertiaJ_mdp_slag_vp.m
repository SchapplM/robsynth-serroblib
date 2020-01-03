% Calculate joint inertia matrix for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5RPRRP7_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:37
% EndTime: 2019-12-31 18:45:38
% DurationCPUTime: 0.26s
% Computational Cost: add. (272->104), mult. (491->155), div. (0->0), fcn. (382->6), ass. (0->39)
t68 = sin(pkin(8));
t61 = pkin(1) * t68 + pkin(6);
t70 = sin(qJ(4));
t92 = t61 * t70;
t72 = cos(qJ(4));
t91 = t61 * t72;
t73 = cos(qJ(3));
t90 = t70 * t73;
t69 = cos(pkin(8));
t62 = -pkin(1) * t69 - pkin(2);
t71 = sin(qJ(3));
t57 = -pkin(3) * t73 - pkin(7) * t71 + t62;
t53 = t70 * t57 + t73 * t91;
t64 = t70 ^ 2;
t66 = t72 ^ 2;
t89 = t64 + t66;
t88 = t70 * MDP(14);
t87 = t72 * MDP(13);
t86 = (MDP(18) - MDP(21));
t85 = t89 * pkin(7);
t84 = -MDP(22) * pkin(4) - MDP(19);
t83 = -2 * qJ(5) * MDP(21) - MDP(16);
t82 = -pkin(4) * t72 - qJ(5) * t70;
t81 = -pkin(4) * t70 + qJ(5) * t72;
t50 = -qJ(5) * t73 + t53;
t56 = t72 * t57;
t51 = -t56 + (pkin(4) + t92) * t73;
t80 = t50 * t72 + t51 * t70;
t79 = t72 * MDP(14) - t70 * MDP(15);
t78 = (-t61 * t90 + t56) * MDP(17) - t53 * MDP(18);
t77 = t70 * MDP(19) - t72 * MDP(21);
t76 = (MDP(22) * qJ(5) - t86) * t72 + (-MDP(17) + t84) * t70;
t67 = t73 ^ 2;
t65 = t71 ^ 2;
t63 = t66 * t71;
t60 = pkin(7) * t90;
t59 = -pkin(3) + t82;
t54 = (t61 - t81) * t71;
t1 = [MDP(1) - 0.2e1 * t62 * t73 * MDP(10) + t67 * MDP(16) + (t50 ^ 2 + t51 ^ 2 + t54 ^ 2) * MDP(22) + (t68 ^ 2 + t69 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * (t62 * MDP(11) + (MDP(6) - t79) * t73) * t71 + 0.2e1 * (t51 * MDP(19) - t50 * MDP(21) - t78) * t73 + 0.2e1 * ((-t50 * t70 + t51 * t72) * MDP(20) + t77 * t54) * t71 + (t66 * MDP(12) - 0.2e1 * t70 * t87 + MDP(5) + 0.2e1 * (t70 * MDP(17) + t72 * MDP(18)) * t61) * t65; (-t54 * t73 + t71 * t80) * MDP(22); MDP(4) + (t65 * t89 + t67) * MDP(22); t63 * MDP(13) + t60 * MDP(17) + (-t54 * t72 + t60) * MDP(19) + t80 * MDP(20) - t54 * t70 * MDP(21) + (pkin(7) * t80 + t54 * t59) * MDP(22) + (-t61 * MDP(11) - t88 + MDP(8) + (pkin(7) * t86 - MDP(15)) * t72) * t73 + (MDP(7) - t61 * MDP(10) + t72 * t70 * MDP(12) - t64 * MDP(13) + (-pkin(3) * t70 - t91) * MDP(17) + (-pkin(3) * t72 + t92) * MDP(18) + t77 * t59) * t71; t63 * MDP(20) + (t64 * MDP(20) + MDP(22) * t85 - MDP(11)) * t71 + (-t59 * MDP(22) + MDP(10) + (MDP(17) + MDP(19)) * t72 - t86 * t70) * t73; MDP(9) + t64 * MDP(12) + (pkin(7) ^ 2 * t89 + t59 ^ 2) * MDP(22) + 0.2e1 * MDP(20) * t85 + 0.2e1 * (MDP(17) * pkin(3) - MDP(19) * t59) * t72 + 0.2e1 * (-MDP(18) * pkin(3) - MDP(21) * t59 + t87) * t70; t56 * MDP(19) + t53 * MDP(21) + (-pkin(4) * t51 + qJ(5) * t50) * MDP(22) + ((-0.2e1 * pkin(4) - t92) * MDP(19) + t83) * t73 + (MDP(20) * t82 + t79) * t71 + t78; t76 * t71; t72 * MDP(15) + MDP(20) * t81 + pkin(7) * t76 + t88; 0.2e1 * pkin(4) * MDP(19) + (pkin(4) ^ 2 + (qJ(5) ^ 2)) * MDP(22) - t83; MDP(20) * t71 * t72 + MDP(19) * t73 + MDP(22) * t51; t70 * t71 * MDP(22); (MDP(22) * pkin(7) + MDP(20)) * t70; t84; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
