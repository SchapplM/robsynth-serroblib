% Calculate joint inertia matrix for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP4_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:02
% EndTime: 2019-03-09 03:13:04
% DurationCPUTime: 0.44s
% Computational Cost: add. (425->135), mult. (670->181), div. (0->0), fcn. (523->6), ass. (0->58)
t93 = sin(pkin(9));
t85 = t93 * pkin(1) + pkin(7);
t132 = pkin(4) + t85;
t95 = sin(qJ(5));
t88 = t95 ^ 2;
t97 = cos(qJ(5));
t90 = t97 ^ 2;
t83 = t88 + t90;
t131 = t83 * MDP(26);
t111 = -MDP(26) * pkin(5) - MDP(23);
t115 = MDP(22) - MDP(25);
t130 = (MDP(26) * qJ(6) - t115) * t95 + (MDP(21) - t111) * t97;
t98 = cos(qJ(3));
t129 = 0.2e1 * t98;
t128 = 0.2e1 * MDP(23);
t99 = -pkin(3) - pkin(8);
t96 = sin(qJ(3));
t127 = t96 * pkin(5);
t126 = t96 * t99;
t94 = cos(pkin(9));
t86 = -t94 * pkin(1) - pkin(2);
t102 = -t96 * qJ(4) + t86;
t71 = t98 * t99 + t102;
t75 = t132 * t96;
t68 = t97 * t71 + t95 * t75;
t76 = t132 * t98;
t89 = t96 ^ 2;
t91 = t98 ^ 2;
t125 = t89 + t91;
t124 = pkin(3) * MDP(15);
t123 = t96 * qJ(6);
t122 = MDP(26) * t99;
t113 = t95 * t71 - t97 * t75;
t66 = t113 - t127;
t121 = t66 * MDP(26);
t109 = t95 * pkin(5) - t97 * qJ(6);
t79 = qJ(4) + t109;
t120 = t79 * MDP(26);
t119 = t85 ^ 2 * MDP(15);
t118 = t97 * MDP(26);
t117 = qJ(4) * MDP(15);
t116 = MDP(21) + MDP(23);
t114 = t97 * t95 * MDP(17);
t112 = -MDP(13) + t124;
t110 = t116 * t97;
t80 = pkin(5) * t97 + t95 * qJ(6);
t65 = t123 + t68;
t64 = t65 * t95 - t66 * t97;
t108 = MDP(10) + t112;
t106 = -MDP(11) + MDP(14) + t117;
t104 = -t95 * MDP(18) - t97 * MDP(19);
t103 = -MDP(21) * t113 - t68 * MDP(22);
t101 = -t115 * t95 + t110;
t78 = t83 * t99;
t77 = t83 * t98;
t74 = -t98 * pkin(3) + t102;
t69 = t80 * t98 + t76;
t1 = [MDP(1) - 0.2e1 * t86 * t98 * MDP(10) + (t65 ^ 2 + t66 ^ 2 + t69 ^ 2) * MDP(26) + (MDP(13) * t129 + MDP(15) * t74) * t74 + (t88 * MDP(16) + 0.2e1 * t114 + t119) * t91 + (MDP(20) + MDP(5) + t119) * t89 + (t93 ^ 2 + t94 ^ 2) * MDP(4) * pkin(1) ^ 2 + ((-t65 * t97 - t66 * t95) * MDP(24) + (t97 * MDP(21) - t95 * MDP(22)) * t76 + (t97 * MDP(23) + t95 * MDP(25)) * t69) * t129 + 0.2e1 * t125 * MDP(12) * t85 + 0.2e1 * (t86 * MDP(11) - t74 * MDP(14) + (MDP(6) + t104) * t98 - t66 * MDP(23) + t65 * MDP(25) + t103) * t96; (-t64 * t98 + t69 * t96) * MDP(26); MDP(4) + t125 * MDP(15) + (t83 * t91 + t89) * MDP(26); t69 * t120 - t64 * MDP(24) + (-pkin(3) * MDP(12) + MDP(7)) * t96 + t110 * t126 + (MDP(8) + qJ(4) * MDP(12) + (t88 - t90) * MDP(17)) * t98 + (-t99 * t121 + t96 * MDP(18) + t76 * MDP(22) - t69 * MDP(25) + (qJ(4) * MDP(21) + t79 * MDP(23)) * t98) * t97 + (-t98 * t97 * MDP(16) - t96 * MDP(19) + t76 * MDP(21) + (-t98 * qJ(4) - t126) * MDP(22) + t69 * MDP(23) + (t79 * t98 + t126) * MDP(25) + t65 * t122) * t95 + (t106 * t98 - t108 * t96) * t85; t77 * MDP(24) + (-t122 * t83 + t108) * t98 + (t115 * t97 + t116 * t95 + t106 + t120) * t96; -0.2e1 * t114 - 0.2e1 * t78 * MDP(24) + t90 * MDP(16) + MDP(9) + t99 ^ 2 * t131 + (-0.2e1 * t97 * MDP(25) + t128 * t95 + t120) * t79 + (-0.2e1 * MDP(13) + t124) * pkin(3) + (0.2e1 * t95 * MDP(21) + 0.2e1 * t97 * MDP(22) + 0.2e1 * MDP(14) + t117) * qJ(4); t64 * MDP(26) + (MDP(15) * t85 + MDP(12) + t101) * t96; -t98 * MDP(15) - t77 * MDP(26); -t83 * MDP(24) + t78 * MDP(26) - t112; MDP(15) + t131; t96 * MDP(20) + (-t113 + 0.2e1 * t127) * MDP(23) + (0.2e1 * t123 + t68) * MDP(25) + (-t66 * pkin(5) + t65 * qJ(6)) * MDP(26) + (MDP(24) * t109 + t104) * t98 + t103; -t130 * t98; t97 * MDP(18) - t95 * MDP(19) - t80 * MDP(24) + t130 * t99; t80 * MDP(26) + t101; MDP(20) + pkin(5) * t128 + 0.2e1 * qJ(6) * MDP(25) + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(26); -t95 * t98 * MDP(24) - t96 * MDP(23) + t121; t98 * t118; (MDP(24) - t122) * t97; -t118; t111; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
