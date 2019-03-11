% Calculate joint inertia matrix for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6RPPRPR3_inertiaJ_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:59
% EndTime: 2019-03-09 01:44:59
% DurationCPUTime: 0.26s
% Computational Cost: add. (339->86), mult. (561->131), div. (0->0), fcn. (574->8), ass. (0->42)
t106 = sin(qJ(4));
t81 = sin(pkin(10));
t83 = cos(pkin(10));
t87 = cos(qJ(4));
t65 = -t83 * t106 - t81 * t87;
t85 = sin(qJ(6));
t86 = cos(qJ(6));
t93 = t86 * MDP(22) - t85 * MDP(23);
t112 = t93 * t65;
t67 = -t81 * t106 + t83 * t87;
t111 = (MDP(19) * t86 - MDP(20) * t85) * t67;
t101 = MDP(23) * t86;
t102 = MDP(22) * t85;
t92 = t101 + t102;
t110 = t92 * (pkin(4) * t81 + pkin(8)) - t85 * MDP(19) - t86 * MDP(20);
t109 = -t106 * MDP(13) - t87 * MDP(14);
t63 = t65 ^ 2;
t64 = t67 ^ 2;
t108 = MDP(7) + (t63 + t64) * MDP(16);
t84 = cos(pkin(9));
t76 = -pkin(1) * t84 - pkin(2);
t97 = -pkin(7) + t76;
t90 = t106 * t97;
t62 = -t106 * qJ(5) + t90;
t91 = (-qJ(5) + t97) * t87;
t57 = t62 * t81 - t83 * t91;
t105 = t57 * t67;
t104 = t85 * t86;
t103 = MDP(16) * pkin(4);
t82 = sin(pkin(9));
t69 = t82 * pkin(1) + qJ(3);
t68 = t106 * pkin(4) + t69;
t100 = MDP(18) * t104;
t96 = t65 * t81 - t67 * t83;
t80 = t86 ^ 2;
t79 = t85 ^ 2;
t75 = -pkin(4) * t83 - pkin(5);
t59 = t83 * t62 + t81 * t91;
t56 = -pkin(5) * t65 - pkin(8) * t67 + t68;
t55 = t56 * t85 + t59 * t86;
t54 = t56 * t86 - t59 * t85;
t1 = [MDP(1) + (t69 ^ 2 + t76 ^ 2) * MDP(7) + (t57 ^ 2 + t59 ^ 2 + t68 ^ 2) * MDP(16) + t63 * MDP(21) + (t82 ^ 2 + t84 ^ 2) * MDP(4) * pkin(1) ^ 2 - 0.2e1 * t65 * t111 + (t80 * MDP(17) - 0.2e1 * t100) * t64 + (MDP(8) * t87 - 0.2e1 * t106 * MDP(9)) * t87 + 0.2e1 * t76 * MDP(5) + 0.2e1 * (t59 * t65 + t105) * MDP(15) + 0.2e1 * (t85 * t105 - t54 * t65) * MDP(22) + 0.2e1 * (t86 * t105 + t55 * t65) * MDP(23) + 0.2e1 * (MDP(6) - t109) * t69; (-t57 * t65 + t59 * t67) * MDP(16); MDP(4) + t108; -t63 * t102 + t76 * MDP(7) + MDP(5) + (-t59 * MDP(16) + (-MDP(15) - t101) * t65) * t65 + (-t57 * MDP(16) + (-MDP(15) - t92) * t67) * t67; 0; t108; -t106 * MDP(11) - MDP(14) * t90 + (t97 * MDP(13) + MDP(10)) * t87 - t93 * t57 + t110 * t65 + (MDP(17) * t104 + (-t79 + t80) * MDP(18) + t92 * t75) * t67 + (t96 * MDP(15) + (-t57 * t83 + t59 * t81) * MDP(16)) * pkin(4); t112 + (t65 * t83 + t67 * t81) * t103 + t109; t87 * MDP(13) - t106 * MDP(14) - t96 * t103 + t93 * t67; 0.2e1 * t100 + t79 * MDP(17) + MDP(12) + (t81 ^ 2 + t83 ^ 2) * MDP(16) * pkin(4) ^ 2 - 0.2e1 * t93 * t75; MDP(16) * t68 - t112; 0; 0; 0; MDP(16); -MDP(21) * t65 + t54 * MDP(22) - t55 * MDP(23) + t111; -t92 * t67; t92 * t65; -t110; t93; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
