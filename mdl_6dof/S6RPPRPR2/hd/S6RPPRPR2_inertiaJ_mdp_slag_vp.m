% Calculate joint inertia matrix for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRPR2_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:38
% EndTime: 2019-03-09 01:42:39
% DurationCPUTime: 0.25s
% Computational Cost: add. (362->91), mult. (632->122), div. (0->0), fcn. (644->8), ass. (0->54)
t111 = cos(qJ(4));
t79 = sin(pkin(10));
t81 = cos(pkin(10));
t84 = sin(qJ(4));
t68 = t111 * t79 + t84 * t81;
t115 = t68 ^ 2;
t114 = -2 * MDP(17);
t113 = pkin(4) + pkin(8);
t80 = sin(pkin(9));
t72 = pkin(1) * t80 + qJ(3);
t112 = pkin(7) + t72;
t63 = t112 * t79;
t64 = t112 * t81;
t61 = t111 * t64 - t84 * t63;
t67 = -t111 * t81 + t79 * t84;
t58 = -t67 * pkin(5) + t61;
t110 = t58 * t67;
t83 = sin(qJ(6));
t85 = cos(qJ(6));
t109 = t83 * t85;
t108 = t79 ^ 2 + t81 ^ 2;
t107 = (MDP(19) * pkin(4));
t82 = cos(pkin(9));
t74 = -pkin(1) * t82 - pkin(2);
t106 = MDP(8) * t74;
t105 = t79 * MDP(6);
t104 = t81 * MDP(5);
t103 = t83 * MDP(25);
t102 = t85 * MDP(26);
t101 = qJ(5) * MDP(19);
t100 = MDP(18) - MDP(15);
t99 = 0.2e1 * t68;
t98 = MDP(21) * t109;
t97 = t108 * MDP(8);
t96 = MDP(17) - t107;
t95 = -MDP(14) + t96;
t70 = -pkin(3) * t81 + t74;
t94 = MDP(22) * t83 + MDP(23) * t85;
t93 = MDP(25) * t85 - t83 * MDP(26);
t92 = t102 + t103;
t60 = t111 * t63 + t84 * t64;
t91 = MDP(16) + t93;
t90 = -t92 - t100;
t89 = -qJ(5) * t68 + t70;
t88 = (-MDP(25) * t113 + MDP(22)) * t85 + (MDP(26) * t113 - MDP(23)) * t83;
t78 = t85 ^ 2;
t77 = t83 ^ 2;
t65 = t67 ^ 2;
t59 = pkin(4) * t67 + t89;
t57 = t68 * pkin(5) + t60;
t56 = t113 * t67 + t89;
t55 = t56 * t85 + t57 * t83;
t54 = -t56 * t83 + t57 * t85;
t1 = [MDP(1) + (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) * MDP(19) + (t80 ^ 2 + t82 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t104 + 0.2e1 * t105 + t106) * t74 + (t70 * MDP(15) - t59 * MDP(18)) * t99 + (MDP(9) + MDP(24)) * t115 + (t77 * MDP(20) + 0.2e1 * t98) * t65 + (0.2e1 * t70 * MDP(14) + t59 * t114 + (-MDP(10) + t94) * t99) * t67 + 0.2e1 * (t60 * t68 - t61 * t67) * MDP(16) + 0.2e1 * (-t85 * t110 + t54 * t68) * MDP(25) + 0.2e1 * (t83 * t110 - t55 * t68) * MDP(26) + (0.2e1 * t108 * MDP(7) + t97 * t72) * t72; (t60 * t67 + t61 * t68) * MDP(19); MDP(4) + t97 + (t65 + t115) * MDP(19); MDP(19) * t59 - t104 + t105 + t106 + (MDP(14) - MDP(17)) * t67 + t90 * t68; 0; MDP(8) + MDP(19); (t100 + t101) * t61 + t95 * t60 + t92 * t58 + (-pkin(4) * MDP(16) + MDP(11) + t88) * t68 + (-MDP(12) + MDP(20) * t109 + (-t77 + t78) * MDP(21) - t91 * qJ(5)) * t67; t95 * t67 + (-t90 + t101) * t68; 0; -0.2e1 * t98 + t78 * MDP(20) + MDP(13) + (t114 + t107) * pkin(4) + (0.2e1 * MDP(18) + t101 + 0.2e1 * t102 + 0.2e1 * t103) * qJ(5); MDP(19) * t60 + t91 * t68; t67 * MDP(19); 0; t96; MDP(19); MDP(24) * t68 + t54 * MDP(25) - t55 * MDP(26) + t94 * t67; t93 * t67; -t92; t88; t93; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
