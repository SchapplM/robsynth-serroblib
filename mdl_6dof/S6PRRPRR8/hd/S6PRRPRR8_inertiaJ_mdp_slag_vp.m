% Calculate joint inertia matrix for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRPRR8_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:40
% EndTime: 2019-03-08 22:40:43
% DurationCPUTime: 0.95s
% Computational Cost: add. (678->203), mult. (1656->290), div. (0->0), fcn. (1780->12), ass. (0->94)
t144 = sin(pkin(7));
t154 = cos(qJ(3));
t196 = t144 * t154;
t150 = sin(qJ(3));
t197 = t144 * t150;
t133 = pkin(9) * t197;
t146 = cos(pkin(7));
t203 = pkin(2) * t154;
t173 = -pkin(3) - t203;
t111 = pkin(4) * t197 + t133 + (-pkin(10) + t173) * t146;
t156 = -pkin(3) - pkin(10);
t171 = -qJ(4) * t150 - pkin(2);
t118 = (t154 * t156 + t171) * t144;
t149 = sin(qJ(5));
t153 = cos(qJ(5));
t105 = t111 * t153 - t118 * t149;
t106 = t111 * t149 + t118 * t153;
t126 = t146 * t153 - t149 * t196;
t209 = t126 * MDP(18) + t105 * MDP(21) - t106 * MDP(22);
t208 = 0.2e1 * t146;
t148 = sin(qJ(6));
t152 = cos(qJ(6));
t207 = MDP(28) * t148 + MDP(29) * t152;
t205 = 0.2e1 * MDP(28);
t204 = 0.2e1 * MDP(29);
t202 = pkin(3) * MDP(15);
t103 = -pkin(5) * t197 - t105;
t201 = t103 * t148;
t200 = t103 * t152;
t145 = sin(pkin(6));
t147 = cos(pkin(6));
t151 = sin(qJ(2));
t155 = cos(qJ(2));
t195 = t146 * t155;
t112 = -t147 * t196 + (t150 * t151 - t154 * t195) * t145;
t124 = -t144 * t145 * t155 + t146 * t147;
t108 = -t112 * t153 + t124 * t149;
t199 = t108 * t153;
t125 = t146 * t149 + t153 * t196;
t198 = t125 * t149;
t194 = t148 * t156;
t193 = t152 * t156;
t128 = pkin(2) * t146 * t150 + pkin(9) * t196;
t191 = MDP(15) * qJ(4);
t190 = MDP(17) * t126;
t189 = MDP(21) * t149;
t109 = t112 * t149 + t124 * t153;
t188 = MDP(22) * t109;
t187 = MDP(22) * t149;
t116 = t126 * t152 + t148 * t197;
t186 = MDP(23) * t116;
t185 = MDP(23) * t152;
t184 = MDP(25) * t116;
t115 = t126 * t148 - t152 * t197;
t183 = MDP(26) * t115;
t182 = MDP(27) * t125;
t123 = t146 * t173 + t133;
t179 = t123 * MDP(15);
t177 = t126 * MDP(22);
t176 = t156 * MDP(22);
t175 = -MDP(10) + MDP(13);
t174 = MDP(11) - MDP(14);
t137 = t146 * qJ(4);
t121 = -t137 - t128;
t172 = MDP(24) * t148 * t152;
t170 = MDP(13) - t202;
t169 = t156 * MDP(21) + MDP(18);
t117 = pkin(4) * t196 - t121;
t168 = (t146 * t203 - t133) * MDP(10) - t128 * MDP(11);
t166 = t125 * MDP(21) + t177;
t165 = MDP(25) * t152 - MDP(26) * t148;
t130 = pkin(5) * t149 - pkin(11) * t153 + qJ(4);
t119 = t130 * t152 - t149 * t194;
t120 = t130 * t148 + t149 * t193;
t164 = MDP(28) * t119 - MDP(29) * t120;
t163 = MDP(28) * t152 - MDP(29) * t148;
t161 = -MDP(17) + t165;
t160 = MDP(21) + t163;
t159 = t148 * MDP(25) + t152 * MDP(26) - pkin(11) * t207;
t104 = pkin(11) * t197 + t106;
t107 = pkin(5) * t125 - pkin(11) * t126 + t117;
t100 = t104 * t152 + t107 * t148;
t99 = -t104 * t148 + t107 * t152;
t158 = MDP(28) * t99 - MDP(29) * t100 + t182 - t183 + t184;
t157 = -MDP(19) + t159;
t143 = t153 ^ 2;
t142 = t152 ^ 2;
t140 = t149 ^ 2;
t139 = t148 ^ 2;
t122 = (-pkin(3) * t154 + t171) * t144;
t113 = t147 * t197 - (-t150 * t195 - t151 * t154) * t145;
t102 = t109 * t152 + t113 * t148;
t101 = -t109 * t148 + t113 * t152;
t1 = [MDP(1) + (t112 ^ 2 + t113 ^ 2 + t124 ^ 2) * MDP(15); t124 * t122 * MDP(15) + (t101 * t125 + t108 * t115) * MDP(28) + (-t102 * t125 + t108 * t116) * MDP(29) + (MDP(3) * t155 - MDP(4) * t151) * t145 + (t146 * t175 + t179) * t112 + (-t121 * MDP(15) - t146 * t174 + t166) * t113 + (t113 * t154 * MDP(12) + (MDP(12) * t112 - MDP(21) * t108 - t188) * t150 + (t150 * t174 + t154 * t175) * t124) * t144; (t121 ^ 2 + t122 ^ 2 + t123 ^ 2) * MDP(15) + t126 ^ 2 * MDP(16) + t146 ^ 2 * MDP(9) + MDP(2) + (-0.2e1 * MDP(24) * t115 + t186) * t116 + (-0.2e1 * MDP(19) * t197 + t182 - 0.2e1 * t183 + 0.2e1 * t184 - 0.2e1 * t190) * t125 + (t103 * t115 + t125 * t99) * t205 + (-t100 * t125 + t103 * t116) * t204 + (t123 * MDP(13) - t121 * MDP(14) + t168) * t208 + 0.2e1 * t166 * t117 + ((t150 * MDP(7) + t154 * MDP(8)) * t208 + 0.2e1 * (-t121 * MDP(12) + t122 * MDP(13)) * t154 + ((MDP(20) + MDP(5)) * t150 ^ 2 + 0.2e1 * (t154 * MDP(10) - t150 * MDP(11)) * pkin(2)) * t144 + 0.2e1 * (t123 * MDP(12) - t122 * MDP(14) + MDP(6) * t196 + t209) * t150) * t144; (t101 * t149 + t148 * t199) * MDP(28) + (-t102 * t149 + t152 * t199) * MDP(29) + (-MDP(10) + t170) * t112 + (MDP(22) * t153 - t174 + t189 + t191) * t113; t133 * MDP(13) + (0.2e1 * t137 + t128) * MDP(14) + (-pkin(3) * t123 - qJ(4) * t121) * MDP(15) + qJ(4) * t177 + (MDP(9) + (-0.2e1 * pkin(3) - t203) * MDP(13)) * t146 + (MDP(21) * qJ(4) + t164) * t125 + (MDP(21) * t117 + t158 - t190) * t149 + ((MDP(12) * qJ(4) + MDP(8)) * t154 + (-pkin(3) * MDP(12) + MDP(7) + (-MDP(19) - t176) * t149) * t150) * t144 + (t126 * MDP(16) + t117 * MDP(22) + t116 * t185 + (-t115 * t152 - t116 * t148) * MDP(24) + (-t115 * t156 + t201) * MDP(28) + (-t116 * t156 + t200) * MDP(29) + t169 * t197 + t161 * t125) * t153 + t168; t140 * MDP(27) + MDP(9) + (-0.2e1 * MDP(13) + t202) * pkin(3) + (0.2e1 * MDP(14) + 0.2e1 * t189 + t191) * qJ(4) + (MDP(23) * t142 + MDP(16) - 0.2e1 * t172) * t143 + (t119 * t149 - t143 * t194) * t205 + (-t120 * t149 - t143 * t193) * t204 + 0.2e1 * (qJ(4) * MDP(22) + t149 * t161) * t153; t112 * MDP(15); t146 * MDP(13) + t179 + (-t115 * t153 - t148 * t198) * MDP(28) + (-t116 * t153 - t152 * t198) * MDP(29) + (MDP(21) * t153 + MDP(12) - t187) * t197; t170 + t207 * (-t140 - t143); MDP(15); -t108 * t160 - t188; MDP(20) * t197 + t148 * t186 + (-t115 * t148 + t116 * t152) * MDP(24) + (-pkin(5) * t115 - t200) * MDP(28) + (-pkin(5) * t116 + t201) * MDP(29) + t157 * t125 + t209; (t157 - t176) * t149 + (t148 * t185 + (-t139 + t142) * MDP(24) + (-pkin(5) * t148 + t193) * MDP(28) + (-pkin(5) * t152 - t194) * MDP(29) + t169) * t153; t153 * t160 - t187; MDP(23) * t139 + 0.2e1 * pkin(5) * t163 + MDP(20) + 0.2e1 * t172; MDP(28) * t101 - MDP(29) * t102; t158; MDP(27) * t149 + t153 * t165 + t164; -t207 * t149; t159; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
