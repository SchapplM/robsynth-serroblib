% Calculate joint inertia matrix for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRPR13_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:29:37
% EndTime: 2019-03-09 11:29:42
% DurationCPUTime: 1.66s
% Computational Cost: add. (1516->269), mult. (3284->390), div. (0->0), fcn. (3464->10), ass. (0->109)
t186 = sin(pkin(6));
t194 = cos(qJ(2));
t231 = t186 * t194;
t191 = sin(qJ(2));
t232 = t186 * t191;
t172 = pkin(8) * t232;
t188 = cos(pkin(6));
t240 = pkin(1) * t194;
t213 = -pkin(2) - t240;
t140 = pkin(3) * t232 + t172 + (-pkin(9) + t213) * t188;
t195 = -pkin(2) - pkin(9);
t212 = -qJ(3) * t191 - pkin(1);
t148 = (t194 * t195 + t212) * t186;
t190 = sin(qJ(4));
t193 = cos(qJ(4));
t135 = t140 * t193 - t190 * t148;
t136 = t140 * t190 + t148 * t193;
t159 = t188 * t193 - t190 * t231;
t251 = t159 * MDP(17) + t135 * MDP(20) - t136 * MDP(21);
t250 = MDP(25) * qJ(5);
t249 = 0.2e1 * t188;
t185 = sin(pkin(11));
t187 = cos(pkin(11));
t229 = t185 ^ 2 + t187 ^ 2;
t248 = t229 * t250;
t246 = 2 * MDP(22);
t245 = 2 * MDP(23);
t244 = 2 * MDP(24);
t243 = -2 * MDP(27);
t242 = 2 * MDP(31);
t241 = 2 * MDP(32);
t239 = pkin(10) * t193;
t238 = pkin(10) + qJ(5);
t237 = pkin(2) * MDP(14);
t236 = pkin(4) * MDP(25);
t130 = -pkin(4) * t232 - t135;
t235 = t130 * t185;
t234 = t130 * t187;
t233 = t185 * t195;
t230 = t190 * t195;
t129 = qJ(5) * t232 + t136;
t161 = t188 * t191 * pkin(1) + pkin(8) * t231;
t178 = t188 * qJ(3);
t151 = -t178 - t161;
t147 = pkin(3) * t231 - t151;
t158 = t188 * t190 + t193 * t231;
t137 = pkin(4) * t158 - qJ(5) * t159 + t147;
t124 = t187 * t129 + t185 * t137;
t167 = pkin(4) * t190 - qJ(5) * t193 + qJ(3);
t150 = t185 * t167 + t187 * t230;
t182 = t190 ^ 2;
t184 = t193 ^ 2;
t228 = -t182 - t184;
t227 = MDP(22) * t187;
t226 = MDP(23) * t185;
t225 = MDP(25) * t190;
t189 = sin(qJ(6));
t192 = cos(qJ(6));
t165 = t185 * t189 - t192 * t187;
t224 = MDP(31) * t165;
t223 = t130 * MDP(25);
t143 = t159 * t185 - t187 * t232;
t144 = t159 * t187 + t185 * t232;
t131 = t192 * t143 + t144 * t189;
t222 = t131 * MDP(29);
t132 = -t143 * t189 + t144 * t192;
t221 = t132 * MDP(26);
t220 = t132 * MDP(28);
t157 = t165 * t193;
t219 = t157 * MDP(26);
t218 = t158 * MDP(30);
t216 = t190 * MDP(30);
t215 = t195 * MDP(21);
t214 = t195 * MDP(25);
t123 = -t129 * t185 + t187 * t137;
t210 = t229 * MDP(24);
t209 = -t123 * t185 + t124 * t187;
t166 = t185 * t192 + t187 * t189;
t207 = -t161 * MDP(10) + (t188 * t240 - t172) * MDP(9);
t205 = -t185 * MDP(22) - t187 * MDP(23);
t155 = t166 * t193;
t204 = -t157 * MDP(28) - t155 * MDP(29);
t154 = t166 * t190;
t156 = t165 * t190;
t203 = -t154 * MDP(31) + t156 * MDP(32);
t202 = -t226 + t227 + t236;
t168 = t238 * t185;
t169 = t238 * t187;
t201 = t166 * MDP(28) - t165 * MDP(29) + (-t168 * t192 - t169 * t189) * MDP(31) - (-t168 * t189 + t169 * t192) * MDP(32);
t200 = MDP(32) * t166 - t202 + t224;
t199 = t143 * MDP(22) + t144 * MDP(23) + t131 * MDP(31) + t132 * MDP(32) + t223;
t198 = (-t143 * t187 + t144 * t185) * MDP(24) + t209 * MDP(25);
t197 = qJ(5) * t205 - MDP(18) + t201;
t177 = -pkin(5) * t187 - pkin(4);
t164 = (pkin(5) * t185 - t195) * t193;
t163 = t187 * t167;
t153 = t188 * t213 + t172;
t152 = (-pkin(2) * t194 + t212) * t186;
t149 = -t185 * t230 + t163;
t141 = -t185 * t239 + t150;
t138 = -t187 * t239 + t163 + (pkin(5) - t233) * t190;
t128 = t138 * t189 + t141 * t192;
t127 = t138 * t192 - t141 * t189;
t125 = pkin(5) * t143 + t130;
t122 = -pkin(10) * t143 + t124;
t121 = pkin(5) * t158 - pkin(10) * t144 + t123;
t120 = t121 * t189 + t122 * t192;
t119 = t121 * t192 - t122 * t189;
t1 = [(t151 ^ 2 + t152 ^ 2 + t153 ^ 2) * MDP(14) + (t123 ^ 2 + t124 ^ 2 + t130 ^ 2) * MDP(25) + t188 ^ 2 * MDP(8) + t159 ^ 2 * MDP(15) + MDP(1) + (t131 * t243 + t221) * t132 + (-0.2e1 * t159 * MDP(16) - 0.2e1 * MDP(18) * t232 + t218 + 0.2e1 * t220 - 0.2e1 * t222) * t158 + (t123 * t158 + t130 * t143) * t246 + (t119 * t158 + t125 * t131) * t242 + (-t124 * t158 + t130 * t144) * t245 + (-t120 * t158 + t125 * t132) * t241 + (-t123 * t144 - t124 * t143) * t244 + (t153 * MDP(12) - t151 * MDP(13) + t207) * t249 + 0.2e1 * (t158 * MDP(20) + t159 * MDP(21)) * t147 + ((t191 * MDP(6) + t194 * MDP(7)) * t249 + 0.2e1 * (-t151 * MDP(11) + t152 * MDP(12)) * t194 + ((MDP(19) + MDP(4)) * t191 ^ 2 + 0.2e1 * (-MDP(10) * t191 + MDP(9) * t194) * pkin(1)) * t186 + 0.2e1 * (t153 * MDP(11) - t152 * MDP(13) + MDP(5) * t231 + t251) * t191) * t186; (-t153 * pkin(2) - t151 * qJ(3)) * MDP(14) + (-t150 * t143 - t149 * t144) * MDP(24) + t172 * MDP(12) + (t123 * t149 + t124 * t150) * MDP(25) + (-t124 * t190 - t150 * t158) * MDP(23) + (t123 * t190 + t149 * t158) * MDP(22) + (qJ(3) * t158 + t147 * t190) * MDP(20) + (t119 * t190 + t125 * t155 + t127 * t158 + t164 * t131) * MDP(31) + (-t120 * t190 - t125 * t157 - t128 * t158 + t164 * t132) * MDP(32) + (-t131 * t190 - t155 * t158) * MDP(29) + (t132 * t190 - t157 * t158) * MDP(28) + (0.2e1 * t178 + t161) * MDP(13) + (t157 * t131 - t132 * t155) * MDP(27) + t158 * t216 - t132 * t219 + ((-0.2e1 * pkin(2) - t240) * MDP(12) + MDP(8)) * t188 + (-t190 * MDP(16) + qJ(3) * MDP(21)) * t159 + (-t158 * MDP(16) + (-t123 * t187 - t124 * t185) * MDP(24) - t130 * t214 + (-t144 * t195 + t234) * MDP(23) + (-t143 * t195 + t235) * MDP(22) + t147 * MDP(21) + t159 * MDP(15)) * t193 + ((qJ(3) * MDP(11) + MDP(7)) * t194 + (-pkin(2) * MDP(11) + MDP(6) + (t195 * MDP(20) + MDP(17)) * t193 + (-MDP(18) - t215) * t190) * t191) * t186 + t207; MDP(8) + t184 * MDP(15) + (t184 * t195 ^ 2 + t149 ^ 2 + t150 ^ 2) * MDP(25) + t182 * MDP(30) - (t155 * t243 - t219) * t157 + (-0.2e1 * MDP(12) + t237) * pkin(2) + (MDP(14) * qJ(3) + 0.2e1 * MDP(21) * t193 + 0.2e1 * MDP(13)) * qJ(3) + 0.2e1 * (-MDP(16) * t193 + qJ(3) * MDP(20) + t204) * t190 + (t149 * t190 - t184 * t233) * t246 + (-t184 * t187 * t195 - t150 * t190) * t245 + (t127 * t190 + t155 * t164) * t242 + (-t128 * t190 - t157 * t164) * t241 + (-t149 * t187 - t150 * t185) * t193 * t244; MDP(11) * t232 + t188 * MDP(12) + t153 * MDP(14) + t203 * t158 + (MDP(20) * t232 - t199) * t193 + (-MDP(21) * t232 + t158 * t205 + t198) * t190; MDP(12) - t237 + t184 * t214 + (-t154 * t190 - t155 * t193) * MDP(31) + (t156 * t190 + t157 * t193) * MDP(32) + (MDP(23) * t228 + t150 * t225) * t187 + (MDP(22) * t228 - t149 * t225) * t185; MDP(14) + (t182 * t229 + t184) * MDP(25); MDP(19) * t232 + (-pkin(4) * t143 - t234) * MDP(22) + (-pkin(4) * t144 + t235) * MDP(23) + t209 * MDP(24) - pkin(4) * t223 + t166 * t221 + (-t131 * t166 - t132 * t165) * MDP(27) + (t125 * t165 + t131 * t177) * MDP(31) + (t125 * t166 + t132 * t177) * MDP(32) + t198 * qJ(5) + t197 * t158 + t251; -t166 * t219 + (-t155 * t166 + t157 * t165) * MDP(27) + (t155 * t177 + t164 * t165) * MDP(31) + (-t157 * t177 + t164 * t166) * MDP(32) + (MDP(17) + t205 * pkin(4) + (MDP(20) + t202) * t195) * t193 + (t197 - t215) * t190 + (MDP(24) + t250) * (-t149 * t185 + t150 * t187); (-MDP(21) + t210 + t248) * t190 + (MDP(20) - t200) * t193; 0.2e1 * t177 * t224 + MDP(19) + (-0.2e1 * t226 + 0.2e1 * t227 + t236) * pkin(4) + (MDP(26) * t166 + t165 * t243 + t177 * t241) * t166 + (0.2e1 * t210 + t248) * qJ(5); t199; t155 * MDP(31) - t157 * MDP(32) + (-t205 - t214) * t193; -t193 * MDP(25); t200; MDP(25); t119 * MDP(31) - t120 * MDP(32) + t218 + t220 - t222; t127 * MDP(31) - t128 * MDP(32) + t204 + t216; t203; t201; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
