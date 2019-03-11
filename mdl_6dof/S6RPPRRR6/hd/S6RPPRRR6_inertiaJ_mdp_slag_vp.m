% Calculate joint inertia matrix for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR6_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:31:40
% EndTime: 2019-03-09 02:31:41
% DurationCPUTime: 0.37s
% Computational Cost: add. (318->116), mult. (569->158), div. (0->0), fcn. (531->6), ass. (0->59)
t109 = cos(qJ(4));
t104 = sin(qJ(6));
t105 = sin(qJ(5));
t107 = cos(qJ(6));
t108 = cos(qJ(5));
t90 = t104 * t108 + t105 * t107;
t80 = t90 * t109;
t89 = t104 * t105 - t107 * t108;
t82 = t89 * t109;
t143 = -t82 * MDP(26) - t80 * MDP(27);
t142 = t109 * MDP(16) + MDP(8);
t135 = pkin(8) + pkin(9);
t92 = t135 * t105;
t93 = t135 * t108;
t117 = t90 * MDP(26) - t89 * MDP(27) + (-t104 * t93 - t107 * t92) * MDP(29) - (-t104 * t92 + t107 * t93) * MDP(30);
t139 = MDP(22) * t105 + MDP(23) * t108;
t141 = t105 * MDP(19) + t108 * MDP(20) - pkin(8) * t139 + t117;
t115 = MDP(22) * t108 - MDP(23) * t105;
t83 = t89 * MDP(29);
t131 = MDP(30) * t90 + t83;
t140 = -t115 + t131;
t137 = -2 * MDP(25);
t136 = 0.2e1 * MDP(30);
t134 = pkin(9) * t109;
t106 = sin(qJ(4));
t133 = t106 * pkin(5);
t102 = pkin(1) + qJ(3);
t79 = t90 * t106;
t81 = t89 * t106;
t132 = -MDP(29) * t79 + MDP(30) * t81;
t101 = -pkin(7) + qJ(2);
t127 = t101 * t108;
t119 = t106 * t127;
t91 = pkin(4) * t106 - pkin(8) * t109 + t102;
t68 = t119 + (t91 - t134) * t105;
t130 = t107 * t68;
t128 = t101 * t105;
t126 = t105 * t108;
t125 = t90 * MDP(24);
t121 = MDP(21) + MDP(28);
t120 = MDP(28) * t106 + t143;
t118 = MDP(18) * t126;
t87 = t108 * t91;
t67 = -t108 * t134 + t87 + (pkin(5) - t128) * t106;
t64 = -t104 * t68 + t107 * t67;
t116 = MDP(19) * t108 - MDP(20) * t105;
t113 = (MDP(29) * t107 - MDP(30) * t104) * pkin(5);
t112 = MDP(15) - t140;
t110 = (qJ(2) ^ 2);
t100 = t109 ^ 2;
t99 = t108 ^ 2;
t98 = t106 ^ 2;
t97 = t105 ^ 2;
t95 = -pkin(5) * t108 - pkin(4);
t88 = (pkin(5) * t105 - t101) * t109;
t74 = t105 * t91 + t119;
t73 = -t106 * t128 + t87;
t65 = t104 * t67 + t130;
t1 = [-(2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2 + t110) * MDP(6)) + (t102 ^ 2 + t110) * MDP(9) + MDP(1) + t121 * t98 - (-MDP(24) * t82 + t137 * t80) * t82 + (MDP(17) * t99 + MDP(10) - 0.2e1 * t118) * t100 + 0.2e1 * (-t100 * t127 - t106 * t74) * MDP(23) + 0.2e1 * (-t100 * t128 + t106 * t73) * MDP(22) + (-t106 * t65 - t82 * t88) * t136 + 0.2e1 * (t106 * t64 + t80 * t88) * MDP(29) + (2 * (MDP(5) + MDP(7)) * qJ(2)) + 0.2e1 * (MDP(15) * t106 + t142) * t102 + 0.2e1 * (t143 + (-MDP(11) + t116) * t109) * t106; -(pkin(1) * MDP(6)) - t102 * MDP(9) - t106 * t112 + MDP(4) - t142; MDP(6) + MDP(9); MDP(7) + qJ(2) * MDP(9) + (-t106 * t79 - t109 * t80) * MDP(29) + (t106 * t81 + t109 * t82) * MDP(30) + t139 * (-t100 - t98); 0; MDP(9); -t82 * t125 + (-t80 * t90 + t82 * t89) * MDP(25) + (t80 * t95 + t88 * t89) * MDP(29) + (-t82 * t95 + t88 * t90) * MDP(30) + (-t101 * MDP(16) - MDP(13) + t141) * t106 + (MDP(12) + t101 * MDP(15) + MDP(17) * t126 + (-t97 + t99) * MDP(18) + (-pkin(4) * t105 + t127) * MDP(22) + (-pkin(4) * t108 - t128) * MDP(23)) * t109; 0; -t106 * MDP(16) + t109 * t112; 0.2e1 * t118 + 0.2e1 * t95 * t83 + t97 * MDP(17) + MDP(14) + 0.2e1 * t115 * pkin(4) + (t136 * t95 + t137 * t89 + t125) * t90; t106 * MDP(21) + t73 * MDP(22) - t74 * MDP(23) + (t107 * t133 + t64) * MDP(29) + (-t130 + (-t67 - t133) * t104) * MDP(30) + t116 * t109 + t120; t140; -t106 * t139 + t132; t141; 0.2e1 * t113 + t121; t64 * MDP(29) - MDP(30) * t65 + t120; t131; t132; t117; MDP(28) + t113; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
