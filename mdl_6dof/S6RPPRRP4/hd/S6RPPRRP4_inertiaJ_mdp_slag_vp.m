% Calculate joint inertia matrix for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRRP4_inertiaJ_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:20
% EndTime: 2019-03-09 02:06:22
% DurationCPUTime: 0.46s
% Computational Cost: add. (496->140), mult. (757->200), div. (0->0), fcn. (656->6), ass. (0->60)
t141 = MDP(27) * pkin(8) + MDP(25);
t102 = sin(qJ(5));
t104 = cos(qJ(5));
t125 = MDP(23) - MDP(26);
t126 = MDP(22) + MDP(24);
t117 = pkin(5) * t104 + qJ(6) * t102;
t89 = -pkin(4) - t117;
t131 = t89 * MDP(27);
t140 = t125 * t102 - t126 * t104 - MDP(15) + t131;
t139 = 0.2e1 * pkin(5);
t138 = pkin(5) * t102;
t105 = cos(qJ(4));
t132 = t104 * t105;
t103 = sin(qJ(4));
t100 = sin(pkin(9));
t101 = cos(pkin(9));
t106 = -pkin(1) - pkin(2);
t86 = qJ(2) * t100 - t101 * t106;
t84 = pkin(3) + t86;
t80 = pkin(4) * t105 + pkin(8) * t103 + t84;
t88 = t101 * qJ(2) + t100 * t106;
t85 = -pkin(7) + t88;
t75 = t102 * t80 + t85 * t132;
t96 = t102 ^ 2;
t98 = t104 ^ 2;
t137 = t96 + t98;
t135 = t102 * t85;
t134 = qJ(6) * t104;
t133 = t102 * t105;
t130 = MDP(18) * t104;
t129 = MDP(25) * t103;
t128 = MDP(27) * t103;
t127 = t105 * MDP(15);
t124 = t137 * pkin(8);
t123 = -MDP(27) * pkin(5) - MDP(24);
t121 = t126 * t102;
t120 = t125 * t104;
t119 = 0.2e1 * qJ(6) * MDP(26) + MDP(21);
t118 = -MDP(22) + t123;
t72 = qJ(6) * t105 + t75;
t79 = t104 * t80;
t73 = -t79 + (-pkin(5) + t135) * t105;
t116 = t102 * t73 + t104 * t72;
t82 = t100 * t133 + t101 * t104;
t83 = t100 * t132 - t101 * t102;
t115 = t102 * t82 + t104 * t83;
t114 = MDP(22) * pkin(4) - MDP(24) * t89;
t113 = -MDP(23) * pkin(4) - MDP(26) * t89;
t112 = MDP(27) * qJ(6) - t125;
t111 = (-t85 * t133 + t79) * MDP(22) - t75 * MDP(23);
t110 = -t104 * MDP(19) + t102 * MDP(20);
t109 = t102 * MDP(19) + t104 * MDP(20);
t108 = t118 * t102;
t99 = t105 ^ 2;
t97 = t103 ^ 2;
t95 = t100 ^ 2;
t93 = t96 * t103;
t90 = t103 * t134;
t76 = t90 + (t85 - t138) * t103;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t86 ^ 2 + t88 ^ 2) * MDP(9) + 0.2e1 * t84 * t127 + t99 * MDP(21) + (t72 ^ 2 + t73 ^ 2 + t76 ^ 2) * MDP(27) + 0.2e1 * t86 * MDP(7) + 0.2e1 * t88 * MDP(8) + 0.2e1 * (-MDP(24) * t73 + MDP(26) * t72 + t111) * t105 + (MDP(17) * t98 - 0.2e1 * t102 * t130 + MDP(10) + 0.2e1 * (-t102 * MDP(22) - t104 * MDP(23)) * t85) * t97 + 0.2e1 * (-t84 * MDP(16) + (MDP(11) + t110) * t105 + (t102 * t72 - t104 * t73) * MDP(25) + (-t102 * MDP(24) + t104 * MDP(26)) * t76) * t103; -MDP(4) - pkin(1) * MDP(6) + (t72 * t83 + t73 * t82) * MDP(27) - t125 * (t100 * t104 * t97 + t105 * t83) - t126 * t82 * t105 + (t102 * t83 - t104 * t82) * t129 + (t88 * MDP(9) - t97 * t121 + t76 * t128 + MDP(8)) * t100 + (t103 * MDP(16) - MDP(9) * t86 - MDP(7) - t127) * t101; MDP(6) + (t101 ^ 2 + t95) * MDP(9) + (t82 ^ 2 + t83 ^ 2 + t97 * t95) * MDP(27); (t116 * t103 - t105 * t76) * MDP(27); (-t100 * t105 + t115) * t128; MDP(9) + (t137 * t97 + t99) * MDP(27); t93 * MDP(18) + (-t104 * MDP(24) - t102 * MDP(26) + t131) * t76 + (-t85 * MDP(16) - MDP(13) + (-t120 - t121) * pkin(8) + t109) * t105 + (-t85 * MDP(15) - t98 * MDP(18) - MDP(12) + (-MDP(22) * t85 - t113) * t104 + (-MDP(17) * t104 + MDP(23) * t85 + t114) * t102) * t103 + t141 * t116; (-t105 * MDP(16) + t140 * t103) * t100 + t141 * t115; t93 * MDP(25) + (t98 * MDP(25) + MDP(27) * t124 - MDP(16)) * t103 - t140 * t105; MDP(14) + t96 * MDP(17) + (t137 * pkin(8) ^ 2 + t89 ^ 2) * MDP(27) + 0.2e1 * MDP(25) * t124 + 0.2e1 * t114 * t104 + 0.2e1 * (t113 + t130) * t102; t79 * MDP(24) + t75 * MDP(26) + (-pkin(5) * t73 + qJ(6) * t72) * MDP(27) + ((t139 - t135) * MDP(24) + t119) * t105 + (t117 * MDP(25) + t110) * t103 + t111; t112 * t83 + t118 * t82; t90 * MDP(27) + (-t120 + t108) * t103; (t134 - t138) * MDP(25) + (t112 * t104 + t108) * pkin(8) + t109; MDP(24) * t139 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(27) + t119; -MDP(24) * t105 + MDP(27) * t73 - t104 * t129; t82 * MDP(27); t102 * t128; t141 * t102; t123; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
