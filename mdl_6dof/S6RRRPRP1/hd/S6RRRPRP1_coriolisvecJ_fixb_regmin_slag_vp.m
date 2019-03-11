% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:33:01
% EndTime: 2019-03-09 16:33:10
% DurationCPUTime: 3.20s
% Computational Cost: add. (8093->367), mult. (20863->486), div. (0->0), fcn. (15482->8), ass. (0->221)
t182 = sin(pkin(10));
t183 = cos(pkin(10));
t188 = cos(qJ(3));
t185 = sin(qJ(3));
t259 = t182 * t185;
t283 = pkin(2) * qJD(3);
t189 = cos(qJ(2));
t295 = pkin(7) + pkin(8);
t160 = t295 * t189;
t155 = qJD(1) * t160;
t143 = t188 * t155;
t186 = sin(qJ(2));
t159 = t295 * t186;
t153 = qJD(1) * t159;
t248 = qJD(1) * t189;
t235 = t188 * t248;
t249 = qJD(1) * t186;
t237 = t185 * t249;
t136 = -t235 + t237;
t265 = qJ(4) * t136;
t98 = t153 * t185 - t143 + t265;
t138 = -t185 * t248 - t188 * t249;
t133 = t138 * qJ(4);
t139 = t185 * t155;
t251 = -t188 * t153 - t139;
t99 = t133 + t251;
t266 = -t182 * t98 - t183 * t99 + (t183 * t188 - t259) * t283;
t242 = qJD(1) * qJD(2);
t234 = t189 * t242;
t241 = qJD(2) + qJD(3);
t112 = qJD(3) * t235 + t188 * t234 - t237 * t241;
t152 = t185 * t189 + t186 * t188;
t298 = qJD(1) * t152;
t194 = t241 * t298;
t192 = t183 * t112 - t182 * t194;
t307 = qJD(5) * t241 + t192;
t184 = sin(qJ(5));
t245 = qJD(5) * t184;
t230 = -t183 * t136 + t138 * t182;
t263 = t230 * t184;
t306 = t245 - t263;
t187 = cos(qJ(5));
t297 = qJD(5) - t230;
t224 = t297 * t187;
t77 = t182 * t112 + t183 * t194;
t305 = -t184 * t77 - t297 * t224;
t304 = -0.2e1 * t242;
t258 = t183 * t185;
t267 = -t182 * t99 + t183 * t98 + (t182 * t188 + t258) * t283;
t282 = qJD(2) * pkin(2);
t145 = -t153 + t282;
t209 = -t145 * t185 - t143;
t238 = qJD(2) * t295;
t220 = qJD(1) * t238;
t147 = t189 * t220;
t128 = t188 * t147;
t146 = t186 * t220;
t228 = t185 * t146 - t128;
t195 = -t112 * qJ(4) + qJD(3) * t209 + t138 * qJD(4) + t228;
t247 = qJD(3) * t185;
t227 = -t185 * t147 - t155 * t247;
t223 = qJD(3) * t145 - t146;
t301 = t188 * t223;
t42 = -qJ(4) * t194 - t136 * qJD(4) + t227 + t301;
t20 = t182 * t42 - t183 * t195;
t229 = t188 * t145 - t139;
t91 = t133 + t229;
t81 = pkin(3) * t241 + t91;
t92 = -t209 - t265;
t86 = t182 * t92;
t47 = t183 * t81 - t86;
t45 = -pkin(4) * t241 - t47;
t303 = t20 * t187 - t45 * t245;
t300 = qJ(6) * t263 + t187 * qJD(6);
t176 = pkin(2) * t249;
t210 = -t136 * t182 - t183 * t138;
t294 = pkin(3) * t138;
t66 = pkin(4) * t210 - pkin(9) * t230 - t294;
t63 = t176 + t66;
t299 = t184 * t63 - t266 * t187;
t96 = t184 * t241 + t187 * t210;
t296 = t96 ^ 2;
t293 = pkin(3) * t182;
t292 = pkin(3) * t183;
t291 = pkin(5) * t187;
t278 = t183 * t92;
t48 = t182 * t81 + t278;
t46 = pkin(9) * t241 + t48;
t175 = -pkin(2) * t189 - pkin(1);
t158 = t175 * qJD(1);
t118 = t136 * pkin(3) + qJD(4) + t158;
t55 = -pkin(4) * t230 - pkin(9) * t210 + t118;
t25 = -t184 * t46 + t187 * t55;
t15 = -qJ(6) * t96 + t25;
t10 = pkin(5) * t297 + t15;
t290 = t10 - t15;
t244 = qJD(5) * t187;
t40 = t184 * t307 + t210 * t244;
t94 = t184 * t210 - t187 * t241;
t289 = -t184 * t40 - t94 * t244;
t53 = t183 * t91 - t86;
t288 = t184 * t66 + t187 * t53;
t105 = -qJ(4) * t152 - t159 * t188 - t160 * t185;
t151 = t185 * t186 - t188 * t189;
t208 = t159 * t185 - t160 * t188;
t106 = -qJ(4) * t151 - t208;
t73 = t105 * t182 + t106 * t183;
t70 = t187 * t73;
t114 = t183 * t151 + t152 * t182;
t115 = -t151 * t182 + t152 * t183;
t207 = pkin(3) * t151 + t175;
t71 = pkin(4) * t114 - pkin(9) * t115 + t207;
t286 = t184 * t71 + t70;
t174 = pkin(2) * t188 + pkin(3);
t132 = pkin(2) * t258 + t182 * t174;
t126 = pkin(9) + t132;
t253 = -qJ(6) - t126;
t226 = qJD(5) * t253;
t285 = t184 * t226 - t299 + t300;
t179 = t187 * qJ(6);
t218 = pkin(5) * t210 - t179 * t230;
t58 = t187 * t63;
t284 = t187 * t226 - t218 - t58 + (-qJD(6) - t266) * t184;
t281 = t10 * t187;
t280 = t210 * t94;
t279 = t230 * t45;
t276 = t184 * t96;
t74 = t187 * t77;
t116 = t241 * t151;
t199 = t152 * qJD(3);
t117 = qJD(2) * t152 + t199;
t79 = -t116 * t183 - t117 * t182;
t275 = t187 * t79;
t274 = t187 * t94;
t273 = t187 * t96;
t39 = -t187 * t307 + t210 * t245;
t271 = t39 * t184;
t270 = t96 * t210;
t171 = pkin(9) + t293;
t252 = -qJ(6) - t171;
t225 = qJD(5) * t252;
t269 = t184 * t225 - t288 + t300;
t232 = t184 * t53 - t187 * t66;
t268 = -t184 * qJD(6) + t187 * t225 - t218 + t232;
t264 = t297 * t210;
t262 = t115 * t184;
t261 = t138 * t136;
t260 = t158 * t138;
t191 = qJD(1) ^ 2;
t256 = t189 * t191;
t190 = qJD(2) ^ 2;
t255 = t190 * t186;
t254 = t190 * t189;
t250 = t186 ^ 2 - t189 ^ 2;
t246 = qJD(3) * t188;
t154 = t186 * t238;
t156 = t189 * t238;
t200 = -t188 * t154 - t185 * t156 - t159 * t246 - t160 * t247;
t59 = -qJ(4) * t117 - qJD(4) * t151 + t200;
t197 = qJD(3) * t208 + t185 * t154 - t188 * t156;
t60 = t116 * qJ(4) - t152 * qJD(4) + t197;
t30 = t182 * t60 + t183 * t59;
t177 = t186 * t282;
t233 = pkin(3) * t117 + t177;
t78 = -t116 * t182 + t183 * t117;
t36 = pkin(4) * t78 - pkin(9) * t79 + t233;
t240 = t184 * t36 + t187 * t30 + t71 * t244;
t172 = -pkin(4) - t292;
t236 = t115 * t244;
t29 = t182 * t59 - t183 * t60;
t52 = t182 * t91 + t278;
t231 = pkin(1) * t304;
t72 = -t183 * t105 + t106 * t182;
t131 = -pkin(2) * t259 + t174 * t183;
t26 = t184 * t55 + t187 * t46;
t221 = t20 * t184 + t210 * t26 + t45 * t244;
t125 = -pkin(4) - t131;
t219 = t306 * pkin(5);
t217 = t20 * t115 - t73 * t77;
t16 = -qJ(6) * t94 + t26;
t216 = -t16 * t184 - t281;
t215 = t210 * t48 + t230 * t47;
t214 = -t126 * t77 - t279;
t213 = -t171 * t77 - t279;
t212 = t274 + t276;
t211 = -qJ(6) * t79 - qJD(6) * t115;
t7 = pkin(5) * t40 + t20;
t206 = -t297 * t306 + t74;
t205 = -t210 * t25 - t303;
t204 = -t187 * t39 - t245 * t96;
t203 = t184 * t79 + t236;
t202 = t158 * t136 - t227;
t21 = t182 * t195 + t183 * t42;
t193 = pkin(3) * t194 + qJD(2) * t176;
t33 = t77 * pkin(4) - pkin(9) * t192 + t193;
t201 = t184 * t33 + t187 * t21 + t55 * t244 - t245 * t46;
t32 = t187 * t33;
t198 = -qJD(5) * t26 - t184 * t21 + t32;
t1 = t77 * pkin(5) + t39 * qJ(6) - t96 * qJD(6) + t198;
t3 = -qJ(6) * t40 - qJD(6) * t94 + t201;
t196 = qJD(5) * t216 - t1 * t184 + t16 * t263 + t3 * t187 + t230 * t281;
t149 = t171 * t187 + t179;
t148 = t252 * t184;
t120 = t126 * t187 + t179;
t119 = t253 * t184;
t100 = -t136 ^ 2 + t138 ^ 2;
t93 = t94 ^ 2;
t85 = -t138 * t241 - t194;
t84 = t136 * t241 + t112;
t69 = t187 * t71;
t37 = t94 * pkin(5) + qJD(6) + t45;
t35 = t187 * t36;
t27 = -qJ(6) * t262 + t286;
t24 = pkin(5) * t114 - t115 * t179 - t184 * t73 + t69;
t17 = t224 * t96 - t271;
t14 = -t270 - t305;
t13 = t206 + t280;
t6 = t212 * t230 + t204 + t289;
t5 = -qJ(6) * t236 + (-qJD(5) * t73 + t211) * t184 + t240;
t4 = t78 * pkin(5) - t184 * t30 + t35 + t211 * t187 + (-t70 + (qJ(6) * t115 - t71) * t184) * qJD(5);
t2 = [0, 0, 0, 0.2e1 * t186 * t234, t250 * t304, t254, -t255, 0, -pkin(7) * t254 + t186 * t231, pkin(7) * t255 + t189 * t231, t112 * t152 + t116 * t138, -t112 * t151 + t116 * t136 + t138 * t117 - t152 * t194, -t116 * t241, -t117 * t241, 0, t136 * t177 + t158 * t117 + t197 * t241 + (t175 * t199 + (t186 * pkin(2) * t151 + t152 * t175) * qJD(2)) * qJD(1), t175 * t112 - t158 * t116 - t200 * t241 + (-t138 + t298) * t177, -t21 * t114 + t192 * t72 + t210 * t29 + t230 * t30 - t47 * t79 - t48 * t78 + t217, t118 * t233 + t193 * t207 + t20 * t72 + t21 * t73 - t47 * t29 + t48 * t30, t115 * t204 + t273 * t79, -t212 * t79 + (t271 - t187 * t40 + (t184 * t94 - t273) * qJD(5)) * t115, t115 * t74 - t39 * t114 + t96 * t78 + (-t115 * t245 + t275) * t297, -t40 * t114 - t203 * t297 - t262 * t77 - t94 * t78, t114 * t77 + t297 * t78 (-t244 * t73 + t35) * t297 + t69 * t77 + (-t244 * t46 + t32) * t114 + t25 * t78 + t29 * t94 + t72 * t40 + t45 * t236 + ((-qJD(5) * t71 - t30) * t297 + (-qJD(5) * t55 - t21) * t114 + t45 * t79 + t217) * t184 -(-t245 * t73 + t240) * t297 - t286 * t77 - t201 * t114 - t26 * t78 + t29 * t96 - t72 * t39 + t45 * t275 + t303 * t115, t24 * t39 - t27 * t40 - t4 * t96 - t5 * t94 + t216 * t79 + (-t1 * t187 - t3 * t184 + (t10 * t184 - t16 * t187) * qJD(5)) * t115, t3 * t27 + t16 * t5 + t1 * t24 + t10 * t4 + t7 * (pkin(5) * t262 + t72) + t37 * (pkin(5) * t203 + t29); 0, 0, 0, -t186 * t256, t250 * t191, 0, 0, 0, t191 * pkin(1) * t186, pkin(1) * t256, -t261, t100, t84, t85, 0, -t136 * t176 - t155 * t246 - t223 * t185 - t128 + t260 + (t143 + (-t153 - t283) * t185) * t241, t138 * t176 + t251 * t241 + (-t241 * t283 - t223) * t188 + t202, -t131 * t192 - t132 * t77 + t267 * t210 + t266 * t230 + t215, t21 * t132 - t20 * t131 - t118 * (t176 - t294) + t266 * t48 - t267 * t47, t17, t6, t14, t13, -t264, t125 * t40 + t267 * t94 + t214 * t184 + (-t126 * t244 - t184 * t266 - t58) * t297 + t205, -t125 * t39 + t267 * t96 + t214 * t187 + (t126 * t245 + t299) * t297 + t221, t119 * t39 - t120 * t40 - t284 * t96 - t285 * t94 + t196, t3 * t120 + t1 * t119 + t7 * (t125 - t291) + (t219 + t267) * t37 + t285 * t16 + t284 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t261, t100, t84, t85, 0, -qJD(2) * t209 + t228 + t260, t229 * t241 + t202 - t301, -t192 * t292 - t210 * t52 - t230 * t53 - t77 * t293 + t215, t47 * t52 - t48 * t53 + (t118 * t138 + t182 * t21 - t183 * t20) * pkin(3), t17, t6, t14, t13, -t264, t172 * t40 - t52 * t94 + t213 * t184 + (-t171 * t244 + t232) * t297 + t205, -t172 * t39 - t52 * t96 + t213 * t187 + (t171 * t245 + t288) * t297 + t221, t148 * t39 - t149 * t40 - t268 * t96 - t269 * t94 + t196, t3 * t149 + t1 * t148 + t7 * (t172 - t291) + (t219 - t52) * t37 + t269 * t16 + t268 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210 ^ 2 - t230 ^ 2, t210 * t47 - t230 * t48 + t193, 0, 0, 0, 0, 0, t206 - t280, -t270 + t305 (t274 - t276) * t230 - t204 + t289, -t37 * t210 + (t16 * t297 + t1) * t187 + (-t10 * t297 + t3) * t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96 * t94, -t93 + t296, t297 * t94 - t39, t297 * t96 - t40, t77, t26 * t297 - t45 * t96 + t198, t25 * t297 + t45 * t94 - t201, pkin(5) * t39 - t290 * t94, t290 * t16 + (-t37 * t96 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93 - t296, t10 * t96 + t16 * t94 + t7;];
tauc_reg  = t2;
