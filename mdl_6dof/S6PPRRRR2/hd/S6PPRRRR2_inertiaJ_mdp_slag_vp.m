% Calculate joint inertia matrix for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR2_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:05:25
% EndTime: 2019-03-08 19:05:27
% DurationCPUTime: 0.51s
% Computational Cost: add. (480->135), mult. (1196->216), div. (0->0), fcn. (1387->14), ass. (0->78)
t133 = sin(qJ(5));
t137 = cos(qJ(5));
t132 = sin(qJ(6));
t136 = cos(qJ(6));
t113 = t132 * t133 - t136 * t137;
t114 = t132 * t137 + t136 * t133;
t169 = pkin(10) + pkin(11);
t117 = t169 * t133;
t118 = t169 * t137;
t148 = t114 * MDP(22) - t113 * MDP(23) + (-t136 * t117 - t132 * t118) * MDP(25) - (-t132 * t117 + t136 * t118) * MDP(26);
t174 = t148 - (MDP(18) * t133 + MDP(19) * t137) * pkin(10) + t133 * MDP(15) + t137 * MDP(16);
t141 = (MDP(25) * t136 - MDP(26) * t132) * pkin(5);
t134 = sin(qJ(4));
t105 = t114 * t134;
t106 = t113 * t134;
t155 = -t106 * MDP(22) - t105 * MDP(23);
t138 = cos(qJ(4));
t116 = -t138 * pkin(4) - t134 * pkin(10) - pkin(3);
t112 = t137 * t116;
t158 = t134 * t137;
t168 = pkin(9) * t133;
t91 = -pkin(11) * t158 + t112 + (-pkin(5) - t168) * t138;
t166 = pkin(9) * t138;
t150 = t137 * t166;
t96 = t150 + (-pkin(11) * t134 + t116) * t133;
t84 = -t132 * t96 + t136 * t91;
t85 = t132 * t91 + t136 * t96;
t173 = t84 * MDP(25) - t85 * MDP(26) + t155;
t171 = -2 * MDP(21);
t170 = 0.2e1 * MDP(26);
t167 = pkin(9) * t137;
t127 = sin(pkin(7));
t128 = sin(pkin(6));
t129 = cos(pkin(13));
t130 = cos(pkin(7));
t131 = cos(pkin(6));
t107 = -t128 * t129 * t127 + t131 * t130;
t126 = sin(pkin(13));
t135 = sin(qJ(3));
t139 = cos(qJ(3));
t161 = t129 * t130;
t163 = t127 * t135;
t93 = t131 * t163 + (t126 * t139 + t135 * t161) * t128;
t89 = t107 * t134 + t93 * t138;
t162 = t127 * t139;
t92 = -t131 * t162 + (t126 * t135 - t139 * t161) * t128;
t80 = -t89 * t133 + t92 * t137;
t81 = t92 * t133 + t89 * t137;
t78 = -t132 * t81 + t136 * t80;
t79 = t132 * t80 + t136 * t81;
t165 = t78 * MDP(25) - t79 * MDP(26);
t109 = t134 * t130 + t138 * t163;
t94 = -t133 * t109 - t137 * t162;
t95 = t137 * t109 - t133 * t162;
t86 = -t132 * t95 + t136 * t94;
t87 = t132 * t94 + t136 * t95;
t164 = t86 * MDP(25) - t87 * MDP(26);
t160 = t133 * t134;
t159 = t133 * t137;
t154 = t113 * MDP(25);
t153 = t114 * MDP(20);
t152 = t134 * MDP(12);
t151 = MDP(17) + MDP(24);
t149 = MDP(14) * t159;
t147 = t137 * MDP(15) - t133 * MDP(16);
t145 = t137 * MDP(18) - t133 * MDP(19);
t142 = t138 * MDP(11) + MDP(4) - t152;
t140 = t114 * MDP(26) - MDP(11) - t145 + t154;
t124 = t137 ^ 2;
t123 = t134 ^ 2;
t122 = t133 ^ 2;
t120 = -t137 * pkin(5) - pkin(4);
t115 = (pkin(5) * t133 + pkin(9)) * t134;
t108 = -t138 * t130 + t134 * t163;
t102 = t133 * t116 + t150;
t101 = -t133 * t166 + t112;
t88 = -t107 * t138 + t93 * t134;
t1 = [MDP(1) + (t131 ^ 2 + (t126 ^ 2 + t129 ^ 2) * t128 ^ 2) * MDP(2); t131 * MDP(2); MDP(2); -t93 * MDP(5) + (-t80 * t138 + t88 * t160) * MDP(18) + (t81 * t138 + t88 * t158) * MDP(19) + (t88 * t105 - t78 * t138) * MDP(25) + (-t88 * t106 + t79 * t138) * MDP(26) - t142 * t92; (t108 * t160 - t94 * t138) * MDP(18) + (t108 * t158 + t95 * t138) * MDP(19) + (t108 * t105 - t86 * t138) * MDP(25) + (-t108 * t106 + t87 * t138) * MDP(26) + (-t135 * MDP(5) + t142 * t139) * t127; -0.2e1 * pkin(3) * t152 + MDP(3) + t151 * t138 ^ 2 - (-t106 * MDP(20) + t105 * t171) * t106 + (t124 * MDP(13) + MDP(6) - 0.2e1 * t149) * t123 + 0.2e1 * (pkin(3) * MDP(11) + (MDP(7) - t147) * t134 - t155) * t138 + 0.2e1 * (-t101 * t138 + t123 * t168) * MDP(18) + 0.2e1 * (t102 * t138 + t123 * t167) * MDP(19) + 0.2e1 * (t115 * t105 - t84 * t138) * MDP(25) + (-t115 * t106 + t85 * t138) * t170; -t89 * MDP(12) + t140 * t88; -t109 * MDP(12) + t140 * t108; -t106 * t153 + (-t114 * t105 + t106 * t113) * MDP(21) + (t120 * t105 + t115 * t113) * MDP(25) + (-t120 * t106 + t115 * t114) * MDP(26) + (-pkin(9) * MDP(12) + MDP(9) - t174) * t138 + (MDP(8) - pkin(9) * MDP(11) + MDP(13) * t159 + (-t122 + t124) * MDP(14) + (-pkin(4) * t133 - t167) * MDP(18) + (-pkin(4) * t137 + t168) * MDP(19)) * t134; 0.2e1 * t149 + 0.2e1 * t120 * t154 + t122 * MDP(13) + MDP(10) + 0.2e1 * t145 * pkin(4) + (t113 * t171 + t120 * t170 + t153) * t114; t80 * MDP(18) - t81 * MDP(19) + t165; t94 * MDP(18) - t95 * MDP(19) + t164; t101 * MDP(18) - t102 * MDP(19) + (-t151 - t141) * t138 + t147 * t134 + t173; t174; 0.2e1 * t141 + t151; t165; t164; -t138 * MDP(24) + t173; t148; MDP(24) + t141; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
