% Calculate joint inertia matrix for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RRPRRP5_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:05:28
% EndTime: 2019-03-09 12:05:31
% DurationCPUTime: 0.97s
% Computational Cost: add. (1594->223), mult. (3761->332), div. (0->0), fcn. (4099->10), ass. (0->102)
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t186 = t166 * MDP(26);
t216 = t163 * MDP(22) + t166 * MDP(23) - (t163 * MDP(25) + t186) * pkin(10) - MDP(16);
t215 = 0.2e1 * MDP(25);
t214 = 0.2e1 * MDP(26);
t213 = 2 * MDP(27);
t168 = cos(qJ(2));
t212 = pkin(1) * t168;
t211 = pkin(8) + qJ(3);
t210 = -qJ(6) - pkin(10);
t209 = MDP(28) * pkin(5);
t162 = cos(pkin(6));
t148 = t162 * t212;
t160 = sin(pkin(6));
t165 = sin(qJ(2));
t203 = t160 * t165;
t132 = pkin(2) * t162 - t211 * t203 + t148;
t184 = pkin(1) * t162 * t165;
t202 = t160 * t168;
t135 = t211 * t202 + t184;
t159 = sin(pkin(11));
t161 = cos(pkin(11));
t121 = t159 * t132 + t161 * t135;
t118 = pkin(9) * t162 + t121;
t137 = t159 * t203 - t161 * t202;
t138 = (t159 * t168 + t161 * t165) * t160;
t144 = (-pkin(2) * t168 - pkin(1)) * t160;
t124 = pkin(3) * t137 - pkin(9) * t138 + t144;
t164 = sin(qJ(4));
t167 = cos(qJ(4));
t114 = -t164 * t118 + t124 * t167;
t110 = -pkin(4) * t137 - t114;
t208 = t110 * t163;
t150 = pkin(2) * t159 + pkin(9);
t207 = t150 * t163;
t206 = t150 * t166;
t205 = t150 * t167;
t154 = t160 ^ 2;
t204 = t154 * t165;
t201 = t162 * MDP(8);
t200 = t163 * t164;
t199 = t164 * t166;
t131 = t138 * t167 + t162 * t164;
t198 = MDP(14) * t131;
t123 = t131 * t166 + t137 * t163;
t197 = MDP(22) * t123;
t122 = t131 * t163 - t137 * t166;
t196 = MDP(23) * t122;
t130 = t138 * t164 - t162 * t167;
t195 = MDP(24) * t130;
t194 = t123 * MDP(20);
t193 = t131 * MDP(13);
t192 = t150 * MDP(19);
t151 = -pkin(2) * t161 - pkin(3);
t191 = t151 * MDP(19);
t153 = -pkin(5) * t166 - pkin(4);
t190 = t153 * MDP(28);
t189 = t163 * MDP(21);
t188 = t163 * MDP(23);
t187 = t166 * MDP(20);
t185 = t167 * MDP(18);
t183 = t166 * t205;
t182 = t166 * t189;
t181 = -MDP(27) * pkin(5) + MDP(22);
t115 = t118 * t167 + t124 * t164;
t111 = pkin(10) * t137 + t115;
t120 = t132 * t161 - t159 * t135;
t117 = -pkin(3) * t162 - t120;
t113 = pkin(4) * t130 - pkin(10) * t131 + t117;
t107 = -t111 * t163 + t166 * t113;
t180 = -t150 * MDP(18) + MDP(15);
t105 = pkin(5) * t130 - qJ(6) * t123 + t107;
t108 = t111 * t166 + t113 * t163;
t106 = -qJ(6) * t122 + t108;
t179 = -t105 * t163 + t106 * t166;
t143 = -pkin(4) * t167 - pkin(10) * t164 + t151;
t139 = t166 * t143;
t125 = -qJ(6) * t199 + t139 + (-pkin(5) - t207) * t167;
t126 = t183 + (-qJ(6) * t164 + t143) * t163;
t178 = -t125 * t163 + t126 * t166;
t145 = t210 * t163;
t146 = t210 * t166;
t177 = -t145 * t163 - t146 * t166;
t128 = -t163 * t205 + t139;
t129 = t143 * t163 + t183;
t176 = t128 * MDP(25) - t129 * MDP(26);
t175 = MDP(25) * t166 - MDP(26) * t163;
t173 = t166 * MDP(22) - MDP(14) - t188;
t171 = (t165 * MDP(6) + t168 * MDP(7)) * t160;
t170 = MDP(25) * t107 - MDP(26) * t108 + t195 - t196 + t197;
t158 = t167 ^ 2;
t157 = t166 ^ 2;
t156 = t164 ^ 2;
t155 = t163 ^ 2;
t152 = t157 * t164;
t142 = pkin(8) * t202 + t184;
t141 = -pkin(8) * t203 + t148;
t140 = (pkin(5) * t163 + t150) * t164;
t119 = t122 * t199;
t109 = pkin(5) * t122 + t110;
t1 = [(t105 ^ 2 + t106 ^ 2 + t109 ^ 2) * MDP(28) + (t120 ^ 2 + t121 ^ 2 + t144 ^ 2) * MDP(12) + t137 ^ 2 * MDP(17) + MDP(1) + (MDP(4) * t165 + 0.2e1 * MDP(5) * t168) * t204 + (0.2e1 * t137 * MDP(15) + t193) * t131 + (-0.2e1 * t122 * MDP(21) + t194) * t123 + (0.2e1 * t171 + t201) * t162 + (-0.2e1 * MDP(16) * t137 + t195 - 0.2e1 * t196 + 0.2e1 * t197 - 0.2e1 * t198) * t130 + 0.2e1 * (t141 * t162 + t154 * t212) * MDP(9) + 0.2e1 * (-pkin(1) * t204 - t142 * t162) * MDP(10) + 0.2e1 * (t114 * t137 + t117 * t130) * MDP(18) + 0.2e1 * (-t115 * t137 + t117 * t131) * MDP(19) + 0.2e1 * (-t120 * t138 - t121 * t137) * MDP(11) + (t107 * t130 + t110 * t122) * t215 + (-t108 * t130 + t110 * t123) * t214 + (-t105 * t123 - t106 * t122) * t213; t201 + t141 * MDP(9) - t142 * MDP(10) + t131 * t191 - t119 * MDP(21) + (-t126 * t122 - t125 * t123) * MDP(27) + (t105 * t125 + t106 * t126 + t109 * t140) * MDP(28) + t171 + (t151 * MDP(18) + t176) * t130 + ((-t137 * t159 - t138 * t161) * MDP(11) + (t120 * t161 + t121 * t159) * MDP(12)) * pkin(2) + (t198 - MDP(18) * t117 + (MDP(16) - t192) * t137 - t170) * t167 + (t193 + t117 * MDP(19) + (t122 * t150 + t208) * MDP(25) + t110 * t186 + (-t105 * t166 - t106 * t163) * MDP(27) + t180 * t137 + (t150 * MDP(26) + t187 - t189) * t123 + t173 * t130) * t164; MDP(8) - 0.2e1 * t151 * t185 + t158 * MDP(24) + (t125 ^ 2 + t126 ^ 2 + t140 ^ 2) * MDP(28) + (t159 ^ 2 + t161 ^ 2) * MDP(12) * pkin(2) ^ 2 + (MDP(20) * t157 + MDP(13) - 0.2e1 * t182) * t156 + (-t128 * t167 + t156 * t207) * t215 + (t129 * t167 + t156 * t206) * t214 + (0.2e1 * t191 + (-t125 * t166 - t126 * t163) * t213 - 0.2e1 * t173 * t167) * t164; t144 * MDP(12) + (-t122 * t167 - t130 * t200) * MDP(25) + (-t123 * t167 - t130 * t199) * MDP(26) + (t123 * t200 - t119) * MDP(27) + (-t109 * t167 + t179 * t164) * MDP(28) + (-t164 * MDP(19) + t185) * t137; (-t140 * t167 + t178 * t164) * MDP(28); MDP(12) + (t158 + (t155 + t157) * t156) * MDP(28); t131 * MDP(15) + t137 * MDP(17) + t114 * MDP(18) - t115 * MDP(19) + t163 * t194 + (-t122 * t163 + t123 * t166) * MDP(21) + (-pkin(4) * t122 - t110 * t166) * MDP(25) + (-pkin(4) * t123 + t208) * MDP(26) + (t122 * t146 - t123 * t145 + t179) * MDP(27) + (t105 * t145 - t106 * t146 + t109 * t153) * MDP(28) + t216 * t130; t152 * MDP(21) + t178 * MDP(27) + (t125 * t145 - t126 * t146 + t140 * t153) * MDP(28) + (-t192 - t216) * t167 + (t163 * t187 - t155 * MDP(21) + (-pkin(4) * t163 - t206) * MDP(25) + (-pkin(4) * t166 + t207) * MDP(26) + (-t145 * t166 + t146 * t163) * MDP(27) + t180) * t164; t152 * MDP(27) + (t155 * MDP(27) + t177 * MDP(28) - MDP(19)) * t164 + (MDP(18) + t175 - t190) * t167; MDP(17) + t155 * MDP(20) + 0.2e1 * t182 + t177 * t213 + (t145 ^ 2 + t146 ^ 2 + t153 ^ 2) * MDP(28) + 0.2e1 * t175 * pkin(4); (-t123 * MDP(27) + t105 * MDP(28)) * pkin(5) + t170; t125 * t209 - MDP(24) * t167 + (t181 * t166 - t188) * t164 + t176; (-t186 + (-MDP(25) - t209) * t163) * t164; t145 * t209 + (-MDP(26) * pkin(10) + MDP(23)) * t166 + (-MDP(25) * pkin(10) + t181) * t163; MDP(28) * pkin(5) ^ 2 + MDP(24); t109 * MDP(28); t140 * MDP(28); -t167 * MDP(28); t190; 0; MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
