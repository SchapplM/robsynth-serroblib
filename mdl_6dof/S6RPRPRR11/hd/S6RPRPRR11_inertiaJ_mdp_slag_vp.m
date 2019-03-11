% Calculate joint inertia matrix for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR11_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:15:17
% EndTime: 2019-03-09 04:15:21
% DurationCPUTime: 1.26s
% Computational Cost: add. (2965->235), mult. (7850->377), div. (0->0), fcn. (9168->14), ass. (0->111)
t192 = sin(pkin(6));
t239 = cos(pkin(12));
t240 = cos(pkin(7));
t213 = t240 * t239;
t191 = sin(pkin(7));
t241 = cos(pkin(6));
t214 = t241 * t191;
t254 = t192 * t213 + t214;
t190 = sin(pkin(12));
t196 = sin(qJ(3));
t198 = cos(qJ(3));
t158 = t196 * t214 + (t190 * t198 + t196 * t213) * t192;
t215 = t192 * t239;
t165 = t191 * t215 - t241 * t240;
t189 = sin(pkin(13));
t193 = cos(pkin(13));
t146 = t158 * t189 + t165 * t193;
t147 = t158 * t193 - t165 * t189;
t195 = sin(qJ(5));
t245 = cos(qJ(5));
t139 = -t195 * t146 + t245 * t147;
t235 = t190 * t192;
t157 = t196 * t235 - t254 * t198;
t253 = t139 * MDP(20) + t157 * MDP(22);
t252 = (MDP(18) * qJ(4));
t138 = t245 * t146 + t147 * t195;
t251 = t138 * MDP(30);
t244 = pkin(1) * t190;
t170 = qJ(2) * t215 + t241 * t244;
t156 = t254 * pkin(9) + t170;
t218 = pkin(1) * t239;
t180 = t241 * t218;
t159 = t241 * pkin(2) + t180 + (-t240 * pkin(9) - qJ(2)) * t235;
t163 = (-pkin(9) * t190 * t191 - t239 * pkin(2) - pkin(1)) * t192;
t206 = t240 * t159 + t163 * t191;
t140 = -t196 * t156 + t206 * t198;
t249 = 2 * MDP(17);
t248 = 2 * MDP(24);
t247 = 2 * MDP(31);
t246 = 2 * MDP(32);
t243 = pkin(10) + qJ(4);
t242 = pkin(3) * MDP(18);
t194 = sin(qJ(6));
t197 = cos(qJ(6));
t131 = t139 * t197 + t157 * t194;
t238 = t131 * t194;
t174 = t245 * t189 + t195 * t193;
t237 = t174 * t194;
t236 = t174 * t197;
t234 = t191 * t196;
t233 = t191 * t198;
t144 = -t159 * t191 + t240 * t163;
t134 = pkin(3) * t157 - qJ(4) * t158 + t144;
t141 = t198 * t156 + t206 * t196;
t136 = -t165 * qJ(4) + t141;
t128 = t189 * t134 + t193 * t136;
t137 = t165 * pkin(3) - t140;
t231 = MDP(18) * t137;
t230 = MDP(25) * t174;
t130 = t139 * t194 - t157 * t197;
t229 = t130 * MDP(29);
t228 = t131 * MDP(28);
t226 = t157 * MDP(11);
t225 = t157 * MDP(21);
t223 = t158 * MDP(10);
t173 = t189 * t195 - t245 * t193;
t222 = t173 * MDP(30);
t221 = t189 * MDP(16);
t220 = t193 * MDP(15);
t219 = t197 * MDP(26);
t183 = -pkin(4) * t193 - pkin(3);
t217 = MDP(27) * t194 * t197;
t216 = t243 * t189;
t127 = t193 * t134 - t136 * t189;
t212 = -t127 * t189 + t128 * t193;
t210 = t197 * MDP(28) - t194 * MDP(29);
t209 = MDP(31) * t197 - MDP(32) * t194;
t208 = t194 * MDP(31) + t197 * MDP(32);
t125 = pkin(4) * t157 - pkin(10) * t147 + t127;
t126 = -pkin(10) * t146 + t128;
t122 = t245 * t125 - t195 * t126;
t123 = t195 * t125 + t245 * t126;
t205 = -MDP(20) + t210;
t204 = MDP(24) + t209;
t203 = -t220 + t221 + t230 - t242;
t202 = t194 * MDP(28) + t197 * MDP(29) - t208 * pkin(11);
t121 = t157 * pkin(11) + t123;
t129 = t146 * pkin(4) + t137;
t124 = t138 * pkin(5) - t139 * pkin(11) + t129;
t118 = -t121 * t194 + t124 * t197;
t119 = t121 * t197 + t124 * t194;
t201 = t118 * MDP(31) - t119 * MDP(32) + t228 - t229 + t251;
t200 = -MDP(22) + t202;
t188 = t197 ^ 2;
t187 = t194 ^ 2;
t185 = t192 ^ 2;
t176 = t243 * t193;
t169 = -qJ(2) * t235 + t180;
t168 = t189 * t240 + t193 * t234;
t167 = -t189 * t234 + t193 * t240;
t161 = t245 * t176 - t195 * t216;
t160 = t176 * t195 + t245 * t216;
t153 = pkin(5) * t173 - pkin(11) * t174 + t183;
t151 = t195 * t167 + t245 * t168;
t150 = -t245 * t167 + t168 * t195;
t149 = t151 * t197 - t194 * t233;
t148 = -t151 * t194 - t197 * t233;
t143 = t153 * t194 + t161 * t197;
t142 = t153 * t197 - t161 * t194;
t120 = -t157 * pkin(5) - t122;
t1 = [-0.2e1 * (t229 + t253) * t138 + (t118 * t247 - t119 * t246 + t129 * t248 + 0.2e1 * t228 + t251) * t138 + 0.2e1 * (-t169 * t190 + t239 * t170) * MDP(6) * t192 + (t127 ^ 2 + t128 ^ 2 + t137 ^ 2) * MDP(18) + (pkin(1) ^ 2 * t185 + t169 ^ 2 + t170 ^ 2) * MDP(7) + 0.2e1 * (t127 * t157 + t137 * t146) * MDP(15) + 0.2e1 * (-t123 * t157 + t129 * t139) * MDP(25) + 0.2e1 * (-t128 * t157 + t137 * t147) * MDP(16) + 0.2e1 * (t141 * t165 + t144 * t158) * MDP(14) + 0.2e1 * (-t140 * t165 + t144 * t157) * MDP(13) + (t131 * MDP(26) - 0.2e1 * t130 * MDP(27) + t120 * t246) * t131 + t139 ^ 2 * MDP(19) - 0.2e1 * t158 * t157 * MDP(9) + t165 ^ 2 * MDP(12) + MDP(1) - 0.2e1 * t165 * t223 + 0.2e1 * t139 * t225 + 0.2e1 * t165 * t226 + (-t127 * t147 - t128 * t146) * t249 + 0.2e1 * (t169 * t241 + t185 * t218) * MDP(4) + 0.2e1 * (-t170 * t241 - t185 * t244) * MDP(5) + t122 * t157 * t248 + t120 * t130 * t247 + t157 ^ 2 * MDP(23) + t158 ^ 2 * MDP(8); (t240 * t157 - t165 * t233) * MDP(13) + (t240 * t158 + t165 * t234) * MDP(14) + (-t146 * t233 + t157 * t167) * MDP(15) + (-t147 * t233 - t157 * t168) * MDP(16) + (-t146 * t168 - t147 * t167) * MDP(17) + (t127 * t167 + t128 * t168 - t137 * t233) * MDP(18) + (-t138 * t233 - t150 * t157) * MDP(24) + (-t139 * t233 - t151 * t157) * MDP(25) + (t130 * t150 + t138 * t148) * MDP(31) + (t131 * t150 - t138 * t149) * MDP(32) + (-t239 * MDP(4) + t190 * MDP(5) - pkin(1) * MDP(7)) * t192; MDP(7) + (t191 ^ 2 * t198 ^ 2 + t167 ^ 2 + t168 ^ 2) * MDP(18); t223 - t226 - t165 * MDP(12) + t140 * MDP(13) - t141 * MDP(14) + (-pkin(3) * t146 - t137 * t193) * MDP(15) + (-pkin(3) * t147 + t137 * t189) * MDP(16) + t212 * MDP(17) - pkin(3) * t231 + (t183 * t138 - t160 * t157) * MDP(24) + (t183 * t139 - t161 * t157) * MDP(25) + (t160 * t130 + t142 * t138) * MDP(31) + (t160 * t131 - t143 * t138) * MDP(32) + ((-t146 * t193 + t147 * t189) * MDP(17) + t212 * MDP(18) + (-t189 * MDP(15) - t193 * MDP(16)) * t157) * qJ(4) + (t129 * MDP(24) + t201 - t253) * t173 + (t139 * MDP(19) + t225 + t129 * MDP(25) + t131 * t219 + (-t130 * t197 - t238) * MDP(27) + t208 * t120 + t205 * t138) * t174; (t148 * t173 + t150 * t237) * MDP(31) + (-t149 * t173 + t150 * t236) * MDP(32) + (-t196 * MDP(14) + (-MDP(24) * t173 + MDP(13) - t203) * t198) * t191 + (MDP(17) + t252) * (-t167 * t189 + t168 * t193); 0.2e1 * t183 * t230 + MDP(12) + (0.2e1 * t220 - 0.2e1 * t221 + t242) * pkin(3) + (MDP(26) * t188 + MDP(19) - 0.2e1 * t217) * t174 ^ 2 + (0.2e1 * t205 * t174 + t183 * t248 + t222) * t173 + (t142 * t173 + t160 * t237) * t247 + (-t143 * t173 + t160 * t236) * t246 + (t249 + t252) * (t189 ^ 2 + t193 ^ 2) * qJ(4); MDP(15) * t146 + MDP(16) * t147 + t139 * MDP(25) + t204 * t138 + t231; -MDP(18) * t233; t204 * t173 + t203; MDP(18); t139 * MDP(21) + t157 * MDP(23) + t122 * MDP(24) - t123 * MDP(25) + MDP(26) * t238 + (-t130 * t194 + t131 * t197) * MDP(27) + (-pkin(5) * t130 - t120 * t197) * MDP(31) + (-pkin(5) * t131 + t120 * t194) * MDP(32) + t200 * t138; -MDP(25) * t151 - t204 * t150; -t161 * MDP(25) - t204 * t160 + t200 * t173 + (MDP(21) + t194 * t219 + (-t187 + t188) * MDP(27) - t208 * pkin(5)) * t174; 0; MDP(26) * t187 + 0.2e1 * pkin(5) * t209 + MDP(23) + 0.2e1 * t217; t201; MDP(31) * t148 - MDP(32) * t149; t142 * MDP(31) - t143 * MDP(32) + t210 * t174 + t222; t209; t202; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
