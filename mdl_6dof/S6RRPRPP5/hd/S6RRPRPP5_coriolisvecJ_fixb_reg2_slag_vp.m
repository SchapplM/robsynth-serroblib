% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:17
% EndTime: 2019-03-09 10:06:28
% DurationCPUTime: 4.22s
% Computational Cost: add. (3901->438), mult. (8755->524), div. (0->0), fcn. (4908->4), ass. (0->230)
t183 = sin(qJ(2));
t270 = qJD(1) * t183;
t337 = -t270 - qJD(4);
t184 = cos(qJ(4));
t182 = sin(qJ(4));
t268 = qJD(2) * t182;
t185 = cos(qJ(2));
t269 = qJD(1) * t185;
t120 = t184 * t269 + t268;
t286 = t337 * t120;
t255 = qJD(1) * qJD(2);
t241 = t183 * t255;
t71 = qJD(4) * t120 - t182 * t241;
t345 = t71 - t286;
t192 = t71 + t286;
t250 = t184 * t270;
t263 = qJD(4) * t184;
t344 = t250 + t263;
t242 = t182 * t269;
t266 = qJD(2) * t184;
t122 = -t242 + t266;
t216 = t122 * t337;
t72 = qJD(2) * t263 - qJD(4) * t242 - t184 * t241;
t328 = t72 - t216;
t5 = -t182 * t345 + t184 * t328;
t312 = pkin(3) + pkin(7);
t297 = t184 * t71;
t28 = t182 * t216 - t297;
t244 = t120 * t263;
t295 = t72 * t182;
t319 = t120 * t250 + t244 + t295;
t343 = t28 + t319;
t165 = t185 * t255;
t311 = pkin(4) + pkin(5);
t342 = t311 * t165;
t340 = qJ(5) * t165 - qJD(5) * t337;
t339 = t165 * t182 - t337 * t344;
t161 = pkin(7) * t165;
t111 = pkin(3) * t165 + t161;
t264 = qJD(4) * t182;
t163 = pkin(2) * t241;
t225 = pkin(8) * t183 - qJ(3) * t185;
t258 = t183 * qJD(3);
t195 = qJD(2) * t225 - t258;
t62 = qJD(1) * t195 + t163;
t187 = -pkin(2) - pkin(8);
t240 = -qJ(3) * t183 - pkin(1);
t112 = t185 * t187 + t240;
t84 = t112 * qJD(1);
t168 = pkin(7) * t270;
t256 = pkin(3) * t270 + qJD(3) + t168;
t89 = qJD(2) * t187 + t256;
t205 = -t111 * t182 - t184 * t62 - t263 * t89 + t264 * t84;
t39 = -t182 * t84 + t184 * t89;
t197 = -t337 * t39 + t205;
t239 = -t111 * t184 + t182 * t62 + t263 * t84 + t264 * t89;
t40 = t182 * t89 + t184 * t84;
t338 = -t337 * t40 - t239;
t267 = qJD(2) * t183;
t336 = t185 * ((t120 * t184 + t122 * t182) * qJD(4) + t295 + t297) - (t120 * t182 - t122 * t184) * t267;
t181 = t185 ^ 2;
t219 = qJD(1) * t181 + t183 * t337;
t262 = qJD(4) * t185;
t245 = t184 * t262;
t288 = t122 * t185;
t335 = qJD(2) * (t182 * t219 - t288) - t337 * t245 + t71 * t183;
t247 = t182 * t262;
t334 = qJD(2) * (t120 * t185 + t184 * t219) + t337 * t247 + t183 * t72;
t333 = -0.2e1 * t255;
t231 = pkin(4) * t165;
t10 = -t231 + t239;
t154 = t337 * qJ(5);
t35 = -t154 + t40;
t332 = t337 * t35 + t10;
t331 = qJ(6) * t72 + qJD(6) * t120;
t170 = pkin(7) * t269;
t128 = pkin(3) * t269 + t170;
t179 = qJD(2) * qJ(3);
t104 = t179 + t128;
t215 = qJ(5) * t122 - t104;
t32 = -t120 * t311 + qJD(6) + t215;
t326 = (qJD(6) + t32) * t122;
t274 = qJD(5) - t39;
t325 = 0.2e1 * t340;
t130 = t337 * t264;
t151 = t184 * t165;
t324 = t151 + t130;
t156 = t337 ^ 2;
t313 = t122 ^ 2;
t323 = -t156 - t313;
t321 = t184 * qJD(5) - t256;
t314 = t120 ^ 2;
t320 = t313 - t314;
t202 = t120 * t269 - t339;
t213 = qJD(2) * t122 + t339;
t281 = t184 * qJ(5);
t317 = t182 * t311 - t281;
t292 = qJ(5) * t182;
t206 = t184 * t311 + t292;
t276 = qJ(6) + t187;
t235 = qJD(4) * t276;
t257 = t184 * qJD(6);
t285 = t182 * t183;
t173 = pkin(2) * t270;
t97 = qJD(1) * t225 + t173;
t49 = t128 * t184 - t182 * t97;
t310 = (-qJ(6) * t285 - t185 * t311) * qJD(1) - t49 - t182 * t235 + t257;
t259 = t182 * qJD(6);
t50 = t128 * t182 + t184 * t97;
t45 = qJ(5) * t269 + t50;
t309 = qJ(6) * t250 + t184 * t235 + t259 - t45;
t308 = -t206 * t337 - t321;
t283 = t182 * t187;
t307 = -t187 * t244 - t283 * t72;
t227 = pkin(4) * t184 + t292;
t306 = t337 * t227 + t321;
t305 = qJ(5) * t72;
t304 = qJD(2) * pkin(2);
t43 = pkin(4) * t120 - t215;
t303 = t122 * t43;
t178 = qJD(2) * qJD(3);
t133 = pkin(7) * t241 - t178;
t94 = -pkin(3) * t241 - t133;
t15 = pkin(4) * t72 + qJ(5) * t71 - qJD(5) * t122 + t94;
t302 = t15 * t182;
t301 = t15 * t184;
t296 = t187 * t71;
t294 = t94 * t182;
t293 = t94 * t184;
t145 = t312 * t183;
t57 = t112 * t184 + t145 * t182;
t291 = t120 * qJ(5);
t289 = t122 * t120;
t287 = t122 * t187;
t284 = t182 * t185;
t282 = t183 * t184;
t280 = t184 * t185;
t189 = qJD(1) ^ 2;
t279 = t185 * t189;
t188 = qJD(2) ^ 2;
t278 = t188 * t183;
t277 = t188 * t185;
t30 = qJ(6) * t122 + t39;
t275 = qJD(5) - t30;
t146 = t312 * t185;
t180 = t183 ^ 2;
t271 = t180 - t181;
t139 = -pkin(2) * t185 + t240;
t105 = qJD(1) * t139;
t265 = qJD(2) * t185;
t261 = qJD(5) * t182;
t260 = t104 * qJD(2);
t52 = qJ(5) * t183 + t57;
t253 = t337 * t283;
t252 = t337 * t184 * t187;
t251 = t182 * t270;
t249 = t182 * t267;
t248 = t187 * t265;
t243 = t337 * t269;
t56 = -t112 * t182 + t145 * t184;
t238 = pkin(1) * t333;
t237 = qJD(3) - t304;
t236 = -t165 + t289;
t234 = -qJD(2) * t120 + t151;
t150 = t183 * t165;
t83 = -t265 * t337 + t150;
t226 = -pkin(4) * t182 + t281;
t33 = pkin(4) * t337 + t274;
t224 = t182 * t33 + t184 * t35;
t223 = -t182 * t39 + t184 * t40;
t129 = t312 * t265;
t172 = pkin(2) * t267;
t78 = t172 + t195;
t23 = -t112 * t263 + t184 * t129 - t145 * t264 - t182 * t78;
t218 = -0.2e1 * qJD(2) * t105;
t113 = t337 * t251;
t214 = t113 + t130 + t234;
t7 = -pkin(5) * t72 - t15;
t212 = -t7 * t182 - t263 * t32;
t211 = -t7 * t184 + t264 * t32;
t204 = -qJ(3) * t265 - t258;
t102 = t172 + t204;
t80 = qJD(1) * t204 + t163;
t210 = pkin(7) * t188 + qJD(1) * t102 + t80;
t209 = t71 * qJ(6) + t239;
t19 = t311 * t337 + t275;
t6 = -t205 + t340;
t2 = t6 + t331;
t31 = qJ(6) * t120 + t40;
t26 = -t154 + t31;
t194 = t209 - t342;
t3 = -t122 * qJD(6) + t194;
t203 = t2 * t182 - t3 * t184 + t344 * t26 + (t251 + t264) * t19;
t22 = -t112 * t264 + t129 * t182 + t145 * t263 + t184 * t78;
t199 = -t183 * t266 - t247;
t198 = t72 + t216;
t14 = qJ(5) * t265 + qJD(5) * t183 + t22;
t25 = t120 * t199 + t280 * t72;
t137 = t168 + t237;
t141 = -t170 - t179;
t191 = -t133 * t185 + (t137 * t185 + (t141 + t170) * t183) * qJD(2);
t159 = t183 * t279;
t144 = -0.2e1 * t150;
t143 = 0.2e1 * t150;
t140 = t271 * t189;
t138 = qJ(3) - t226;
t136 = t276 * t184;
t135 = t276 * t182;
t134 = t187 * t151;
t127 = t312 * t267;
t125 = -qJ(3) * t269 + t173;
t108 = -qJ(3) - t317;
t103 = t271 * t333;
t88 = t105 * t270;
t77 = t185 * t227 + t146;
t60 = pkin(4) * t122 + t291;
t59 = -t185 * t206 - t146;
t53 = -pkin(4) * t183 - t56;
t47 = -pkin(4) * t269 - t49;
t46 = qJ(6) * t280 + t52;
t44 = -t122 * t311 - t291;
t41 = qJ(6) * t284 - t183 * t311 - t56;
t37 = (t285 * t337 - t288) * qJD(1) + t324;
t36 = (qJD(4) * t226 + t261) * t185 + (-t227 - t312) * t267;
t27 = (qJD(4) * t317 - t261) * t185 + (t206 + t312) * t267;
t24 = t71 * t284 + (-t245 + t249) * t122;
t18 = -pkin(4) * t265 - t23;
t9 = qJ(6) * t199 + t185 * t257 + t14;
t8 = -qJ(6) * t249 + (qJ(6) * t263 - qJD(2) * t311 + t259) * t185 - t23;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t103, t277, t144, -t278, 0, -pkin(7) * t277 + t183 * t238, pkin(7) * t278 + t185 * t238, 0, 0, 0, -t277, t278, t143, t103, t144, t191, t183 * t218 + t185 * t210, -t183 * t210 + t185 * t218, pkin(7) * t191 + t105 * t102 + t80 * t139, t24, t336, -t335, t25, -t334, t83, -t127 * t120 + t146 * t72 - t23 * t337 + (-t184 * t260 - t239) * t183 + (-t104 * t264 + t293 + (qJD(1) * t56 + t39) * qJD(2)) * t185, -t127 * t122 - t146 * t71 + t22 * t337 + (t182 * t260 + t205) * t183 + (-t104 * t263 - t294 + (-qJD(1) * t57 - t40) * qJD(2)) * t185, -t22 * t120 - t23 * t122 + t56 * t71 - t57 * t72 + t223 * t267 + (t205 * t184 - t239 * t182 + (t182 * t40 + t184 * t39) * qJD(4)) * t185, -t104 * t127 + t146 * t94 - t205 * t57 + t22 * t40 + t23 * t39 - t239 * t56, t24, -t335, -t336, t83, t334, t25, t36 * t120 + t18 * t337 + t77 * t72 + (-t266 * t43 - t10) * t183 + (-t43 * t264 + t301 + (-qJD(1) * t53 - t33) * qJD(2)) * t185, -t14 * t120 + t18 * t122 - t52 * t72 - t53 * t71 + t224 * t267 + (-t10 * t182 - t184 * t6 + (t182 * t35 - t184 * t33) * qJD(4)) * t185, -t36 * t122 - t14 * t337 + t77 * t71 + (-t268 * t43 + t6) * t183 + (t43 * t263 + t302 + (qJD(1) * t52 + t35) * qJD(2)) * t185, t10 * t53 + t14 * t35 + t15 * t77 + t18 * t33 + t36 * t43 + t52 * t6, t24, -t336, t335, t25, -t334, t83, -t27 * t120 + t8 * t337 - t59 * t72 + (t266 * t32 - t3) * t183 + ((-qJD(1) * t41 - t19) * qJD(2) + t211) * t185, t27 * t122 - t9 * t337 - t59 * t71 + (t268 * t32 + t2) * t183 + ((qJD(1) * t46 + t26) * qJD(2) + t212) * t185, t9 * t120 - t8 * t122 + t41 * t71 + t46 * t72 + (-t182 * t19 - t184 * t26) * t267 + (t182 * t3 + t184 * t2 + (-t182 * t26 + t184 * t19) * qJD(4)) * t185, t19 * t8 + t2 * t46 + t26 * t9 + t27 * t32 + t3 * t41 + t59 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159, t140, 0, t159, 0, 0, t189 * pkin(1) * t183, pkin(1) * t279, 0, 0, 0, 0, 0, -t159, t140, t159 ((-t141 - t179) * t183 + (-t137 + t237) * t185) * qJD(1), -t125 * t269 + t88, 0.2e1 * t178 + (t105 * t185 + t125 * t183) * qJD(1), -t133 * qJ(3) - t141 * qJD(3) - t105 * t125 + (-t141 * t183 + (-t137 - t304) * t185) * qJD(1) * pkin(7), t28, -t5, t37, t319, t202, t243, qJ(3) * t72 + t49 * t337 + t294 + t134 + t256 * t120 + (t104 * t184 + t253) * qJD(4) + (t104 * t282 - t185 * t39) * qJD(1), -qJ(3) * t71 - t50 * t337 + t293 + t256 * t122 + (-t104 * t182 + t252) * qJD(4) + (t185 * t40 + (-t104 * t183 - t248) * t182) * qJD(1), t50 * t120 + t49 * t122 + (t296 - t338) * t184 + (t39 * t270 + t205 + (t39 + t287) * qJD(4)) * t182 + t307, t94 * qJ(3) - t39 * t49 - t40 * t50 + t256 * t104 + (qJD(4) * t223 - t182 * t205 - t184 * t239) * t187, t28, t37, t5, t243, -t202, t319, t138 * t72 + t302 - t47 * t337 + t134 - t306 * t120 + (t184 * t43 + t253) * qJD(4) + (t185 * t33 + t282 * t43) * qJD(1), t45 * t120 - t47 * t122 + (t296 + t332) * t184 + (-t33 * t270 - t6 + (-t33 + t287) * qJD(4)) * t182 + t307, t138 * t71 - t301 + t45 * t337 + t306 * t122 + (t182 * t43 - t252) * qJD(4) + (-t185 * t35 + (t183 * t43 + t248) * t182) * qJD(1), t15 * t138 - t33 * t47 - t35 * t45 - t306 * t43 + (qJD(4) * t224 - t10 * t184 + t6 * t182) * t187, t28, t5, t122 * t269 - t113 - t324, t319, t202, t243, -t108 * t72 - t310 * t337 + t308 * t120 + (-t32 * t282 + (qJD(2) * t136 + t19) * t185) * qJD(1) + t212, -t108 * t71 - t309 * t337 - t308 * t122 + (-t32 * t285 + (qJD(2) * t135 - t26) * t185) * qJD(1) - t211, t120 * t309 + t122 * t310 + t135 * t72 - t136 * t71 + t203, t7 * t108 + t2 * t135 - t3 * t136 - t19 * t310 + t26 * t309 - t308 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, -t180 * t189 - t188, qJD(2) * t141 + t161 + t88, 0, 0, 0, 0, 0, 0, -t156 * t182 + t234, -t213, -t343, -t182 * t197 + t184 * t338 - t260, 0, 0, 0, 0, 0, 0, t214, -t343, t213, -t43 * qJD(2) - t332 * t184 + (-t33 * t337 + t6) * t182, 0, 0, 0, 0, 0, 0, t214, t213, t343, qJD(2) * t32 + t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, t320, -t192, -t289, -t198, t165, -t104 * t122 + t338, t104 * t120 + t197, 0, 0, t289, -t192, -t320, t165, t198, -t289, -t120 * t60 + 0.2e1 * t231 - t303 + t338, pkin(4) * t71 - t305 + (t35 - t40) * t122 + (t33 - t274) * t120, -t120 * t43 + t122 * t60 - t197 + t325, -t10 * pkin(4) + t6 * qJ(5) + t274 * t35 - t33 * t40 - t43 * t60, t289, -t320, t192, -t289, -t198, t165, t44 * t120 - t31 * t337 - t209 + t326 + 0.2e1 * t342, t120 * t32 - t122 * t44 + t30 * t337 - t205 + t325 + t331, t305 - t311 * t71 + (-t26 + t31) * t122 + (-t19 + t275) * t120, t2 * qJ(5) - t19 * t31 + t26 * t275 - t3 * t311 - t32 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, -t192, t323, t303 + t332, 0, 0, 0, 0, 0, 0, t236, t323, t192, t26 * t337 + t194 - t326; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t328, -t345, -t313 - t314, -t120 * t26 + t122 * t19 + t7;];
tauc_reg  = t1;
