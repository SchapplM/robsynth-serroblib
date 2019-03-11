% Calculate joint inertia matrix for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRP12_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:50:43
% EndTime: 2019-03-09 06:50:48
% DurationCPUTime: 1.65s
% Computational Cost: add. (3430->290), mult. (8994->436), div. (0->0), fcn. (10247->12), ass. (0->113)
t272 = MDP(32) * pkin(11) + MDP(30);
t271 = MDP(27) + MDP(29);
t234 = MDP(28) - MDP(31);
t197 = sin(pkin(12));
t199 = sin(pkin(6));
t200 = cos(pkin(12));
t252 = t199 * t200;
t202 = cos(pkin(6));
t267 = pkin(1) * t202;
t177 = qJ(2) * t252 + t197 * t267;
t198 = sin(pkin(7));
t201 = cos(pkin(7));
t251 = t200 * t201;
t233 = t199 * t251;
t160 = (t198 * t202 + t233) * pkin(9) + t177;
t205 = sin(qJ(3));
t208 = cos(qJ(3));
t189 = t200 * t267;
t255 = t197 * t199;
t163 = pkin(2) * t202 + t189 + (-pkin(9) * t201 - qJ(2)) * t255;
t169 = (-pkin(9) * t197 * t198 - pkin(2) * t200 - pkin(1)) * t199;
t219 = t163 * t201 + t169 * t198;
t150 = -t205 * t160 + t219 * t208;
t270 = 2 * MDP(29);
t269 = 2 * MDP(30);
t193 = t199 ^ 2;
t268 = pkin(1) * t193;
t254 = t198 * t205;
t162 = t202 * t254 + (t197 * t208 + t205 * t251) * t199;
t175 = -t198 * t252 + t201 * t202;
t204 = sin(qJ(4));
t207 = cos(qJ(4));
t155 = t162 * t204 - t175 * t207;
t266 = pkin(5) * t155;
t203 = sin(qJ(5));
t265 = pkin(10) * t203;
t206 = cos(qJ(5));
t264 = pkin(10) * t206;
t262 = pkin(3) * MDP(20);
t261 = qJ(6) * t155;
t154 = -t163 * t198 + t201 * t169;
t253 = t198 * t208;
t161 = -t202 * t253 + t205 * t255 - t208 * t233;
t145 = pkin(3) * t161 - pkin(10) * t162 + t154;
t151 = t160 * t208 + t219 * t205;
t148 = pkin(10) * t175 + t151;
t140 = t145 * t207 - t204 * t148;
t138 = -pkin(4) * t161 - t140;
t260 = t138 * t203;
t259 = t138 * t206;
t156 = t162 * t207 + t175 * t204;
t152 = t156 * t203 - t161 * t206;
t258 = t152 * t206;
t153 = t156 * t206 + t161 * t203;
t257 = t153 * t203;
t178 = -t207 * t201 + t204 * t254;
t256 = t178 * t204;
t250 = t203 * t207;
t141 = t145 * t204 + t148 * t207;
t139 = pkin(11) * t161 + t141;
t147 = -pkin(3) * t175 - t150;
t144 = pkin(4) * t155 - pkin(11) * t156 + t147;
t135 = t206 * t139 + t203 * t144;
t184 = -pkin(4) * t207 - pkin(11) * t204 - pkin(3);
t172 = t203 * t184 + t207 * t264;
t194 = t203 ^ 2;
t196 = t206 ^ 2;
t249 = t194 + t196;
t248 = MDP(30) * t204;
t247 = t152 * MDP(25);
t246 = t153 * MDP(22);
t245 = t153 * MDP(24);
t244 = t155 * MDP(26);
t243 = t156 * MDP(16);
t242 = t156 * MDP(17);
t241 = t161 * MDP(19);
t224 = -pkin(5) * t206 - qJ(6) * t203;
t183 = -pkin(4) + t224;
t240 = t183 * MDP(32);
t239 = t203 * MDP(24);
t238 = t204 * MDP(21);
t237 = t206 * MDP(22);
t236 = t206 * MDP(23);
t232 = -pkin(10) * MDP(20) + MDP(17);
t231 = -pkin(10) * MDP(21) + MDP(18);
t230 = -MDP(32) * pkin(5) - MDP(29);
t229 = t139 * t203 - t206 * t144;
t226 = -0.2e1 * qJ(6) * MDP(31) - MDP(26);
t225 = -MDP(27) + t230;
t223 = -pkin(5) * t203 + qJ(6) * t206;
t222 = MDP(32) * qJ(6) - t234;
t132 = t135 + t261;
t133 = t229 - t266;
t221 = t132 * t206 + t133 * t203;
t167 = -qJ(6) * t207 + t172;
t181 = t206 * t184;
t168 = -t181 + (pkin(5) + t265) * t207;
t217 = t167 * t206 + t168 * t203;
t216 = t206 * MDP(24) - t203 * MDP(25);
t215 = t206 * MDP(25) + t239;
t214 = (-pkin(10) * t250 + t181) * MDP(27) - t172 * MDP(28);
t213 = t203 * MDP(29) - t206 * MDP(31);
t212 = -MDP(16) + t216;
t211 = t168 * MDP(29) - t167 * MDP(31) - t214;
t210 = -MDP(27) * t229 - t135 * MDP(28) + t244 + t245 - t247;
t190 = pkin(11) * t250;
t179 = t201 * t204 + t207 * t254;
t176 = -qJ(2) * t255 + t189;
t174 = (pkin(10) - t223) * t204;
t165 = t179 * t206 - t203 * t253;
t164 = t179 * t203 + t206 * t253;
t136 = pkin(5) * t152 - qJ(6) * t153 + t138;
t1 = [t175 ^ 2 * MDP(12) + (pkin(1) ^ 2 * t193 + t176 ^ 2 + t177 ^ 2) * MDP(7) + t156 ^ 2 * MDP(15) + (t132 ^ 2 + t133 ^ 2 + t136 ^ 2) * MDP(32) + MDP(1) + (0.2e1 * t175 * MDP(10) + MDP(8) * t162) * t162 + (-0.2e1 * t152 * MDP(23) + t246) * t153 + (-0.2e1 * t175 * MDP(11) - 0.2e1 * t162 * MDP(9) + t241 + 0.2e1 * t242) * t161 + (-0.2e1 * t161 * MDP(18) - 0.2e1 * t243 + t244 + 0.2e1 * t245 - 0.2e1 * t247) * t155 + 0.2e1 * (t150 * t175 + t154 * t161) * MDP(13) + 0.2e1 * (-t151 * t175 + t154 * t162) * MDP(14) + 0.2e1 * (-t141 * t161 + t147 * t156) * MDP(21) + 0.2e1 * (t140 * t161 + t147 * t155) * MDP(20) + 0.2e1 * (t132 * t155 - t136 * t153) * MDP(31) + 0.2e1 * (t138 * t152 - t155 * t229) * MDP(27) + 0.2e1 * (-t135 * t155 + t138 * t153) * MDP(28) + (-t133 * t155 + t136 * t152) * t270 + (-t132 * t152 + t133 * t153) * t269 + 0.2e1 * (-t177 * t202 - t197 * t268) * MDP(5) + 0.2e1 * (t176 * t202 + t200 * t268) * MDP(4) + 0.2e1 * (-t176 * t197 + t177 * t200) * MDP(6) * t199; (t161 * t201 + t175 * t253) * MDP(13) + (t162 * t201 - t175 * t254) * MDP(14) + (-t155 * t253 - t161 * t178) * MDP(20) + (-t156 * t253 - t161 * t179) * MDP(21) + (-t152 * t165 + t153 * t164) * MDP(30) + (t132 * t165 + t133 * t164 + t136 * t178) * MDP(32) + (-MDP(4) * t200 + MDP(5) * t197 - MDP(7) * pkin(1)) * t199 + t271 * (t178 * t152 - t155 * t164) + t234 * (t153 * t178 - t155 * t165); MDP(7) + (t164 ^ 2 + t165 ^ 2 + t178 ^ 2) * MDP(32); t162 * MDP(10) - t161 * MDP(11) + t175 * MDP(12) + t150 * MDP(13) - t151 * MDP(14) - pkin(3) * t156 * MDP(21) + (-t152 * t167 + t153 * t168) * MDP(30) + (t132 * t167 + t133 * t168) * MDP(32) + (t152 * MDP(29) - t153 * MDP(31) + t136 * MDP(32)) * t174 + (-t211 - t262) * t155 + (-t147 * MDP(20) + t133 * MDP(29) - t132 * MDP(31) + t231 * t161 - t210 + t243) * t207 + (t156 * MDP(15) + t147 * MDP(21) + t153 * t237 + (-t257 - t258) * MDP(23) + (pkin(10) * t152 + t260) * MDP(27) + (pkin(10) * t153 + t259) * MDP(28) + (-t132 * t203 + t133 * t206) * MDP(30) + t232 * t161 + t213 * t136 + t212 * t155) * t204; (t164 * t168 + t165 * t167 + t174 * t178) * MDP(32) + (t164 * t206 - t165 * t203) * t248 + t234 * (t165 * t207 + t206 * t256) + (-t205 * MDP(14) + (t207 * MDP(20) + MDP(13) - t238) * t208) * t198 + t271 * (t164 * t207 + t203 * t256); MDP(12) - 0.2e1 * pkin(3) * t238 + (t167 ^ 2 + t168 ^ 2 + t174 ^ 2) * MDP(32) + 0.2e1 * ((-t167 * t203 + t168 * t206) * MDP(30) + t213 * t174) * t204 + (t207 * MDP(26) - 0.2e1 * t212 * t204 + 0.2e1 * t211 + 0.2e1 * t262) * t207 + (t196 * MDP(22) - 0.2e1 * t203 * t236 + MDP(15) + 0.2e1 * (t203 * MDP(27) + t206 * MDP(28)) * pkin(10)) * t204 ^ 2; t242 + t241 + t140 * MDP(20) - t141 * MDP(21) + t203 * t246 + (-t152 * t203 + t153 * t206) * MDP(23) + (-pkin(4) * t152 - t259) * MDP(27) + (-pkin(4) * t153 + t260) * MDP(28) + (-t136 * t206 + t152 * t183) * MDP(29) + t221 * MDP(30) + (-t136 * t203 - t153 * t183) * MDP(31) + t136 * t240 + ((t257 - t258) * MDP(30) + t221 * MDP(32)) * pkin(11) + (-MDP(18) + (-t203 * t271 - t234 * t206) * pkin(11) + t215) * t155; -t179 * MDP(21) + (t234 * t203 - t206 * t271 - MDP(20) + t240) * t178 + t272 * (t164 * t203 + t165 * t206); t190 * MDP(27) + (-t174 * t206 + t190) * MDP(29) + t217 * MDP(30) - t174 * t203 * MDP(31) + (t217 * pkin(11) + t174 * t183) * MDP(32) + (-t239 + (t234 * pkin(11) - MDP(25)) * t206 + t231) * t207 + (t203 * t237 + (-t194 + t196) * MDP(23) + (-pkin(4) * t203 - t264) * MDP(27) + (-pkin(4) * t206 + t265) * MDP(28) + t213 * t183 + t232) * t204; MDP(19) + t194 * MDP(22) + (t249 * pkin(11) ^ 2 + t183 ^ 2) * MDP(32) + t249 * pkin(11) * t269 + 0.2e1 * (pkin(4) * MDP(27) - t183 * MDP(29)) * t206 + 0.2e1 * (-pkin(4) * MDP(28) - t183 * MDP(31) + t236) * t203; (-t229 + 0.2e1 * t266) * MDP(29) + (-pkin(5) * t153 - qJ(6) * t152) * MDP(30) + (t135 + 0.2e1 * t261) * MDP(31) + (-pkin(5) * t133 + qJ(6) * t132) * MDP(32) + t210; t225 * t164 + t222 * t165; t181 * MDP(29) + t172 * MDP(31) + (-pkin(5) * t168 + qJ(6) * t167) * MDP(32) + ((-0.2e1 * pkin(5) - t265) * MDP(29) + t226) * t207 + (t224 * MDP(30) + t216) * t204 + t214; t223 * MDP(30) + (t225 * t203 + t222 * t206) * pkin(11) + t215; pkin(5) * t270 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32) - t226; -t155 * MDP(29) + t153 * MDP(30) + t133 * MDP(32); t164 * MDP(32); t207 * MDP(29) + MDP(32) * t168 + t206 * t248; t272 * t203; t230; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
