% Calculate joint inertia matrix for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRRPP2_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:52:42
% EndTime: 2019-03-08 22:52:44
% DurationCPUTime: 0.62s
% Computational Cost: add. (574->187), mult. (1132->241), div. (0->0), fcn. (1079->8), ass. (0->66)
t156 = MDP(19) + MDP(23);
t177 = MDP(17) + t156;
t128 = sin(qJ(3));
t176 = 0.2e1 * t128;
t152 = MDP(22) + MDP(26);
t127 = sin(qJ(4));
t133 = pkin(4) + pkin(5);
t130 = cos(qJ(4));
t164 = t130 * qJ(5);
t174 = t133 * t127 - t164;
t173 = pkin(9) - qJ(6);
t172 = MDP(22) * pkin(9);
t170 = cos(pkin(6));
t126 = sin(pkin(6));
t129 = sin(qJ(2));
t169 = t126 * t129;
t132 = cos(qJ(2));
t168 = t126 * t132;
t167 = t127 * t128;
t131 = cos(qJ(3));
t166 = t127 * t131;
t165 = t128 * t130;
t163 = t131 * qJ(5);
t113 = -t131 * pkin(3) - t128 * pkin(9) - pkin(2);
t162 = pkin(8) * t166 - t130 * t113;
t104 = t130 * t131 * pkin(8) + t127 * t113;
t123 = t127 ^ 2;
t125 = t130 ^ 2;
t161 = t123 + t125;
t160 = MDP(11) * t128;
t150 = t127 * qJ(5) + pkin(3);
t109 = t133 * t130 + t150;
t159 = t109 * MDP(26);
t112 = -t130 * pkin(4) - t150;
t158 = t112 * MDP(22);
t157 = t130 * MDP(13);
t155 = MDP(20) - MDP(25);
t154 = MDP(21) - MDP(18);
t153 = MDP(21) + MDP(24);
t122 = t131 * pkin(4);
t101 = t122 + t162;
t151 = -0.2e1 * t163 + t104;
t148 = MDP(18) - t153;
t147 = -pkin(4) * MDP(22) - MDP(19);
t146 = MDP(20) + t172;
t100 = t104 - t163;
t145 = -pkin(4) * t127 + t164;
t143 = -t162 * MDP(17) - t104 * MDP(18);
t142 = -t127 * MDP(23) + t130 * MDP(24);
t140 = -t133 * MDP(26) - MDP(23) + t147;
t114 = t173 * t127;
t115 = t173 * t130;
t139 = -t127 * MDP(14) + t114 * MDP(23) - t115 * MDP(24);
t138 = pkin(3) * MDP(17) - t112 * MDP(19) + t109 * MDP(23) - t115 * MDP(25);
t137 = -pkin(3) * MDP(18) - t112 * MDP(21) + t109 * MDP(24) - t114 * MDP(25);
t135 = qJ(5) ^ 2;
t116 = qJ(6) * t167;
t108 = t170 * t128 + t131 * t169;
t106 = t128 * t169 - t170 * t131;
t105 = (pkin(8) - t145) * t128;
t99 = (-pkin(8) - t174) * t128;
t98 = t108 * t130 - t127 * t168;
t97 = t108 * t127 + t130 * t168;
t92 = t100 + t116;
t91 = t131 * pkin(5) - qJ(6) * t165 + t101;
t1 = [MDP(1) + t152 * (t106 ^ 2 + t97 ^ 2 + t98 ^ 2); (t98 * t100 + t97 * t101 + t106 * t105) * MDP(22) + (-t106 * t99 + t97 * t91 + t98 * t92) * MDP(26) - t155 * (t127 * t98 - t130 * t97) * t128 + (-t129 * MDP(4) + (MDP(10) * t131 + MDP(3) - t160) * t132) * t126 + (-MDP(24) - t154) * (t106 * t165 + t98 * t131) + t177 * (t106 * t167 + t97 * t131); MDP(2) - 0.2e1 * pkin(2) * t160 + (t100 ^ 2 + t101 ^ 2 + t105 ^ 2) * MDP(22) + (t91 ^ 2 + t92 ^ 2 + t99 ^ 2) * MDP(26) + (0.2e1 * pkin(2) * MDP(10) + t131 * MDP(16) + (-t130 * MDP(14) + t127 * MDP(15) + MDP(6)) * t176) * t131 + 0.2e1 * (t101 * MDP(19) - t100 * MDP(21) + t91 * MDP(23) - t92 * MDP(24) - t143) * t131 + ((-t100 * t127 + t101 * t130) * MDP(20) + (t127 * t92 - t130 * t91) * MDP(25) + t142 * t99 + (t127 * MDP(19) - t130 * MDP(21)) * t105) * t176 + (t125 * MDP(12) - 0.2e1 * t127 * t157 + MDP(5) + 0.2e1 * (t127 * MDP(17) + t130 * MDP(18)) * pkin(8)) * t128 ^ 2; -t108 * MDP(11) + (t97 * t114 + t98 * t115) * MDP(26) + (t148 * t127 - t130 * t177 - MDP(10) + t158 - t159) * t106 + (t155 + t172) * (t97 * t127 + t98 * t130); t105 * t158 + (t99 * t109 + t91 * t114 + t92 * t115) * MDP(26) + (MDP(17) + MDP(19)) * pkin(9) * t166 + (-t105 * MDP(19) + t99 * MDP(23) - t92 * MDP(25) + t100 * t146) * t130 + (-t105 * MDP(21) + t99 * MDP(24) - t91 * MDP(25) + t101 * t146) * t127 + (-pkin(8) * MDP(11) + MDP(8) + (-t154 * pkin(9) - MDP(15)) * t130 + t139) * t131 + (MDP(7) - pkin(8) * MDP(10) + (-t123 + t125) * MDP(13) + (-pkin(8) * MDP(17) + t137) * t130 + (t130 * MDP(12) + pkin(8) * MDP(18) - t138) * t127) * t128; MDP(9) + t123 * MDP(12) + (t161 * pkin(9) ^ 2 + t112 ^ 2) * MDP(22) + (t109 ^ 2 + t114 ^ 2 + t115 ^ 2) * MDP(26) + 0.2e1 * t161 * MDP(20) * pkin(9) + 0.2e1 * t138 * t130 + 0.2e1 * (t137 + t157) * t127; (-MDP(17) + t140) * t97 + (t152 * qJ(5) - t148) * t98; (-0.2e1 * t122 - t162) * MDP(19) + t151 * MDP(21) + (-t101 * pkin(4) + t100 * qJ(5)) * MDP(22) - t101 * MDP(23) + (t116 + t151) * MDP(24) + (t92 * qJ(5) - t91 * t133) * MDP(26) + (-MDP(16) + (-pkin(5) - t133) * MDP(23)) * t131 + ((-t155 * qJ(5) - MDP(15)) * t127 + (-pkin(4) * MDP(20) + qJ(6) * MDP(23) + t133 * MDP(25) + MDP(14)) * t130) * t128 + t143; t130 * MDP(15) + t145 * MDP(20) + t174 * MDP(25) + (t115 * qJ(5) - t114 * t133) * MDP(26) + ((qJ(5) * MDP(22) + t154) * t130 + (-MDP(17) + t147) * t127) * pkin(9) - t139; MDP(16) + 0.2e1 * pkin(4) * MDP(19) + (pkin(4) ^ 2 + t135) * MDP(22) + 0.2e1 * t133 * MDP(23) + (t133 ^ 2 + t135) * MDP(26) + 0.2e1 * t153 * qJ(5); t152 * t97; t101 * MDP(22) + t91 * MDP(26) + t156 * t131 + t155 * t165; t114 * MDP(26) + (-MDP(25) + t146) * t127; t140; t152; -t106 * MDP(26); t99 * MDP(26) + t128 * t142; t130 * MDP(23) + t127 * MDP(24) + t159; 0; 0; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
