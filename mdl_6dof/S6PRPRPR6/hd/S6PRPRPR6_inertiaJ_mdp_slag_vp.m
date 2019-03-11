% Calculate joint inertia matrix for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRPR6_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:49:36
% EndTime: 2019-03-08 19:49:37
% DurationCPUTime: 0.47s
% Computational Cost: add. (390->136), mult. (811->208), div. (0->0), fcn. (832->10), ass. (0->71)
t152 = qJ(5) * MDP(18);
t103 = sin(pkin(11));
t105 = cos(pkin(11));
t143 = t103 ^ 2 + t105 ^ 2;
t151 = t143 * t152;
t150 = MDP(17) + t152;
t149 = -2 * MDP(20);
t148 = 2 * MDP(25);
t147 = (pkin(2) * MDP(7));
t111 = cos(qJ(4));
t146 = pkin(9) * t111;
t145 = pkin(9) + qJ(5);
t108 = sin(qJ(4));
t113 = -pkin(2) - pkin(8);
t138 = t108 * t113;
t92 = t108 * pkin(4) - t111 * qJ(5) + qJ(3);
t80 = t103 * t92 + t105 * t138;
t144 = pkin(4) * MDP(18);
t142 = qJ(3) * MDP(7);
t141 = t103 * t113;
t104 = sin(pkin(6));
t109 = sin(qJ(2));
t140 = t104 * t109;
t112 = cos(qJ(2));
t139 = t104 * t112;
t107 = sin(qJ(6));
t110 = cos(qJ(6));
t90 = t107 * t103 - t110 * t105;
t84 = t90 * t111;
t137 = t84 * MDP(19);
t136 = t90 * MDP(24);
t101 = t108 ^ 2;
t102 = t111 ^ 2;
t135 = -t101 - t102;
t134 = MDP(18) * t108;
t132 = t103 * MDP(16);
t131 = t105 * MDP(15);
t130 = t111 * MDP(14);
t129 = t113 * MDP(18);
t128 = MDP(5) - t147;
t127 = t143 * MDP(17);
t106 = cos(pkin(6));
t86 = t106 * t111 - t108 * t139;
t75 = -t86 * t103 + t105 * t140;
t76 = t103 * t140 + t86 * t105;
t125 = -t75 * t103 + t76 * t105;
t91 = t110 * t103 + t107 * t105;
t82 = t91 * t111;
t123 = -t84 * MDP(21) - t82 * MDP(22);
t122 = (-t107 * t76 + t110 * t75) * MDP(24) - (t107 * t75 + t110 * t76) * MDP(25);
t121 = t82 * MDP(24) - t84 * MDP(25);
t120 = -t103 * MDP(15) - t105 * MDP(16);
t119 = t131 - t132 + t144;
t118 = -t120 - t129;
t93 = t145 * t103;
t94 = t145 * t105;
t117 = t91 * MDP(21) - t90 * MDP(22) + (-t107 * t94 - t110 * t93) * MDP(24) - (-t107 * t93 + t110 * t94) * MDP(25);
t116 = t91 * MDP(25) - t119 + t136;
t115 = MDP(13) - t116;
t97 = -t105 * pkin(5) - pkin(4);
t89 = (pkin(5) * t103 - t113) * t111;
t88 = t105 * t92;
t85 = t106 * t108 + t111 * t139;
t83 = t90 * t108;
t81 = t91 * t108;
t79 = -t103 * t138 + t88;
t74 = -t103 * t146 + t80;
t73 = -t105 * t146 + t88 + (pkin(5) - t141) * t108;
t70 = t107 * t73 + t110 * t74;
t69 = -t107 * t74 + t110 * t73;
t1 = [MDP(1) + (t106 ^ 2 + (t109 ^ 2 + t112 ^ 2) * t104 ^ 2) * MDP(7) + (t75 ^ 2 + t76 ^ 2 + t85 ^ 2) * MDP(18); (t75 * t79 + t76 * t80) * MDP(18) + t121 * t85 + (t75 * MDP(15) - t76 * MDP(16) + t122) * t108 + ((-t103 * t76 - t105 * t75) * MDP(17) + t118 * t85) * t111 + ((MDP(3) - t128) * t112 + (t108 * MDP(13) - MDP(4) + MDP(6) + t130 + t142) * t109) * t104; MDP(2) + t102 * MDP(8) + (t102 * t113 ^ 2 + t79 ^ 2 + t80 ^ 2) * MDP(18) + t101 * MDP(23) - (t82 * t149 - t137) * t84 + ((-2 * MDP(5) + t147) * pkin(2)) + (0.2e1 * MDP(6) + 0.2e1 * t130 + t142) * qJ(3) + 0.2e1 * (qJ(3) * MDP(13) - t111 * MDP(9) + t123) * t108 + 0.2e1 * (-t102 * t141 + t79 * t108) * MDP(15) + 0.2e1 * (-t102 * t113 * t105 - t80 * t108) * MDP(16) + 0.2e1 * (t69 * t108 + t89 * t82) * MDP(24) + (-t70 * t108 - t89 * t84) * t148 + 0.2e1 * (-t103 * t80 - t105 * t79) * MDP(17) * t111; -MDP(7) * t139 + (t125 * t108 - t85 * t111) * MDP(18); t102 * t129 + (-t81 * t108 - t111 * t82) * MDP(24) + (t83 * t108 + t111 * t84) * MDP(25) + (t135 * MDP(16) + t80 * t134) * t105 + (t135 * MDP(15) - t79 * t134) * t103 + t128; MDP(7) + (t143 * t101 + t102) * MDP(18); -t86 * MDP(14) - t115 * t85 + t150 * t125; -t91 * t137 + (-t91 * t82 + t84 * t90) * MDP(20) + (t97 * t82 + t89 * t90) * MDP(24) + (-t97 * t84 + t89 * t91) * MDP(25) + (MDP(10) + t120 * pkin(4) + (MDP(13) + t119) * t113) * t111 + (-t113 * MDP(14) + t120 * qJ(5) - MDP(11) + t117) * t108 + t150 * (-t79 * t103 + t80 * t105); (-MDP(14) + t127 + t151) * t108 + t115 * t111; 0.2e1 * t97 * t136 + MDP(12) + (0.2e1 * t131 - 0.2e1 * t132 + t144) * pkin(4) + (MDP(19) * t91 + t97 * t148 + t90 * t149) * t91 + (0.2e1 * t127 + t151) * qJ(5); t85 * MDP(18); t118 * t111 + t121; -t111 * MDP(18); t116; MDP(18); t122; t108 * MDP(23) + t69 * MDP(24) - t70 * MDP(25) + t123; -t81 * MDP(24) + t83 * MDP(25); t117; 0; MDP(23);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
