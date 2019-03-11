% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:10:02
% EndTime: 2019-03-09 10:10:13
% DurationCPUTime: 5.10s
% Computational Cost: add. (9053->378), mult. (23635->529), div. (0->0), fcn. (18567->10), ass. (0->224)
t222 = sin(pkin(10));
t224 = cos(pkin(10));
t227 = sin(qJ(2));
t229 = cos(qJ(2));
t193 = -t222 * t227 + t224 * t229;
t181 = t193 * qJD(1);
t195 = t222 * t229 + t224 * t227;
t183 = t195 * qJD(1);
t226 = sin(qJ(4));
t307 = cos(qJ(4));
t138 = t307 * t181 - t183 * t226;
t134 = qJD(6) - t138;
t221 = sin(pkin(11));
t223 = cos(pkin(11));
t225 = sin(qJ(6));
t228 = cos(qJ(6));
t196 = t221 * t228 + t223 * t225;
t265 = qJD(1) * qJD(2);
t259 = t229 * t265;
t260 = t227 * t265;
t173 = -t222 * t260 + t224 * t259;
t182 = t195 * qJD(2);
t236 = qJD(1) * t182;
t242 = -t226 * t181 - t307 * t183;
t233 = qJD(4) * t242 - t226 * t173 - t307 * t236;
t275 = t228 * t223;
t277 = t221 * t225;
t194 = -t275 + t277;
t343 = t134 * t194;
t344 = -t343 * t134 - t196 * t233;
t342 = t134 * t196;
t326 = t138 * t221;
t341 = pkin(5) * t326;
t340 = pkin(9) * t326;
t218 = qJD(2) + qJD(4);
t286 = t242 * t218;
t339 = t233 - t286;
t254 = -t342 * t134 + t194 * t233;
t122 = t218 * t221 - t223 * t242;
t120 = -t223 * t218 - t221 * t242;
t322 = t228 * t120;
t74 = t122 * t225 + t322;
t297 = t242 * t74;
t338 = t254 - t297;
t325 = t138 * t223;
t302 = -qJ(3) - pkin(7);
t207 = t302 * t229;
t201 = qJD(1) * t207;
t188 = t222 * t201;
t206 = t302 * t227;
t200 = qJD(1) * t206;
t298 = qJD(2) * pkin(2);
t192 = t200 + t298;
t140 = t224 * t192 + t188;
t304 = pkin(8) * t183;
t116 = qJD(2) * pkin(3) + t140 - t304;
t276 = t224 * t201;
t141 = t222 * t192 - t276;
t305 = pkin(8) * t181;
t119 = t141 + t305;
t67 = t226 * t116 + t307 * t119;
t64 = qJ(5) * t218 + t67;
t264 = -pkin(2) * t229 - pkin(1);
t251 = t264 * qJD(1);
t203 = qJD(3) + t251;
t146 = -t181 * pkin(3) + t203;
t71 = -pkin(4) * t138 + qJ(5) * t242 + t146;
t38 = -t221 * t64 + t223 * t71;
t258 = qJD(2) * t302;
t178 = t229 * qJD(3) + t227 * t258;
t162 = t178 * qJD(1);
t179 = -t227 * qJD(3) + t229 * t258;
t163 = t179 * qJD(1);
t124 = -t162 * t222 + t224 * t163;
t104 = -pkin(8) * t173 + t124;
t125 = t224 * t162 + t222 * t163;
t105 = -pkin(8) * t236 + t125;
t261 = qJD(4) * t307;
t269 = qJD(4) * t226;
t232 = t226 * t104 + t307 * t105 + t116 * t261 - t119 * t269;
t27 = t218 * qJD(5) + t232;
t211 = pkin(2) * t260;
t147 = pkin(3) * t236 + t211;
t95 = t307 * t173 + t181 * t261 - t183 * t269 - t226 * t236;
t37 = -pkin(4) * t233 - t95 * qJ(5) + qJD(5) * t242 + t147;
t9 = t221 * t37 + t223 * t27;
t7 = t9 * t223;
t337 = t38 * t325 + t7;
t336 = -pkin(5) * t242 - pkin(9) * t325;
t29 = -t307 * t104 + t226 * t105 + t116 * t269 + t119 * t261;
t296 = t221 * t95;
t15 = pkin(5) * t296 + t29;
t16 = -pkin(5) * t138 - pkin(9) * t122 + t38;
t39 = t221 * t71 + t223 * t64;
t24 = -pkin(9) * t120 + t39;
t4 = t16 * t225 + t228 * t24;
t66 = t307 * t116 - t226 * t119;
t63 = -t218 * pkin(4) + qJD(5) - t66;
t55 = t120 * pkin(5) + t63;
t335 = t15 * t196 - t4 * t242 - t343 * t55;
t3 = t16 * t228 - t225 * t24;
t334 = t15 * t194 + t3 * t242 + t342 * t55;
t247 = -qJD(6) * t122 - t296;
t266 = qJD(6) * t228;
t288 = -t120 * t266 + t95 * t275;
t22 = t225 * t247 + t288;
t246 = t120 * t225 - t122 * t228;
t333 = t22 * t196 + t343 * t246;
t293 = t246 * t242;
t332 = -t293 + t344;
t23 = -qJD(6) * t246 + t196 * t95;
t331 = -t22 * t194 - t196 * t23 + t342 * t246 + t343 * t74;
t283 = t138 * t218;
t329 = t95 - t283;
t328 = t134 * t242;
t327 = t134 * t246;
t321 = t242 * t138;
t320 = -t138 ^ 2 + t242 ^ 2;
t8 = -t221 * t27 + t223 * t37;
t319 = t138 * t39 - t8;
t318 = t29 * t221 - t242 * t39;
t315 = t146 * t242 - t29;
t314 = -t223 * t29 + t242 * t38;
t313 = -t146 * t138 - t232;
t97 = -pkin(4) * t242 - qJ(5) * t138;
t312 = -0.2e1 * t265;
t144 = -t200 * t222 + t276;
t126 = t144 - t305;
t145 = t224 * t200 + t188;
t127 = t145 - t304;
t213 = pkin(2) * t224 + pkin(3);
t306 = pkin(2) * t222;
t238 = t226 * t213 + t307 * t306;
t289 = t238 * qJD(4) + t307 * t126 - t226 * t127;
t209 = t226 * t306;
t245 = -qJD(4) * t209 + t213 * t261;
t161 = qJD(5) + t245;
t78 = t226 * t126 + t307 * t127;
t270 = qJD(1) * t227;
t155 = pkin(2) * t270 + pkin(3) * t183;
t79 = t155 + t97;
t43 = t221 * t79 + t223 * t78;
t311 = -t161 * t223 + t43;
t50 = t221 * t97 + t223 * t66;
t310 = -qJD(5) * t223 + t50;
t150 = t224 * t206 + t207 * t222;
t131 = -pkin(8) * t195 + t150;
t151 = t222 * t206 - t224 * t207;
t132 = pkin(8) * t193 + t151;
t309 = t307 * t131 - t226 * t132;
t308 = -qJD(6) + t134;
t303 = t223 * pkin(5);
t217 = t223 * pkin(9);
t129 = -t178 * t222 + t224 * t179;
t185 = t193 * qJD(2);
t107 = -pkin(8) * t185 + t129;
t130 = t224 * t178 + t222 * t179;
t108 = -pkin(8) * t182 + t130;
t44 = qJD(4) * t309 + t226 * t107 + t307 * t108;
t241 = t307 * t193 - t226 * t195;
t100 = t241 * qJD(4) - t226 * t182 + t307 * t185;
t143 = t226 * t193 + t307 * t195;
t101 = t143 * qJD(4) + t307 * t182 + t226 * t185;
t216 = t227 * t298;
t156 = pkin(3) * t182 + t216;
t48 = pkin(4) * t101 - qJ(5) * t100 - qJD(5) * t143 + t156;
t12 = t221 * t48 + t223 * t44;
t165 = -pkin(3) * t193 + t264;
t90 = -pkin(4) * t241 - qJ(5) * t143 + t165;
t92 = t226 * t131 + t307 * t132;
t54 = t221 * t90 + t223 * t92;
t294 = t223 * t95;
t290 = t289 - t341;
t287 = t100 * t221;
t280 = t143 * t221;
t279 = t143 * t223;
t231 = qJD(1) ^ 2;
t274 = t229 * t231;
t230 = qJD(2) ^ 2;
t273 = t230 * t227;
t272 = t230 * t229;
t271 = t227 ^ 2 - t229 ^ 2;
t267 = qJD(6) * t143;
t2 = -pkin(5) * t233 - pkin(9) * t294 + t8;
t5 = -pkin(9) * t296 + t9;
t263 = t228 * t2 - t225 * t5;
t11 = -t221 * t44 + t223 * t48;
t49 = -t221 * t66 + t223 * t97;
t42 = -t221 * t78 + t223 * t79;
t53 = -t221 * t92 + t223 * t90;
t256 = pkin(1) * t312;
t253 = t225 * t2 + t228 * t5;
t252 = -t221 * t8 + t7;
t250 = t221 * t38 - t223 * t39;
t31 = -pkin(5) * t241 - pkin(9) * t279 + t53;
t40 = -pkin(9) * t280 + t54;
t249 = -t225 * t40 + t228 * t31;
t248 = t225 * t31 + t228 * t40;
t177 = -t307 * t213 - pkin(4) + t209;
t176 = qJ(5) + t238;
t148 = (-pkin(9) - t176) * t221;
t244 = -qJD(6) * t148 + t311 - t340;
t149 = t176 * t223 + t217;
t243 = qJD(6) * t149 + t161 * t221 + t336 + t42;
t204 = (-pkin(9) - qJ(5)) * t221;
t240 = -qJD(6) * t204 + t310 - t340;
t205 = qJ(5) * t223 + t217;
t239 = qJD(5) * t221 + qJD(6) * t205 + t336 + t49;
t237 = t100 * t63 + t143 * t29 - t309 * t95;
t235 = -pkin(4) * t95 + qJ(5) * t233 - (-qJD(5) + t63) * t138;
t234 = t176 * t233 + t177 * t95 + (t161 - t63) * t138;
t45 = t92 * qJD(4) - t307 * t107 + t226 * t108;
t214 = -pkin(4) - t303;
t157 = t177 - t303;
t99 = t194 * t143;
t98 = t196 * t143;
t58 = pkin(5) * t280 - t309;
t56 = t67 + t341;
t36 = t100 * t196 + t266 * t279 - t267 * t277;
t35 = -t100 * t194 - t196 * t267;
t21 = pkin(5) * t287 + t45;
t10 = -pkin(9) * t287 + t12;
t6 = pkin(5) * t101 - t100 * t217 + t11;
t1 = [0, 0, 0, 0.2e1 * t227 * t259, t271 * t312, t272, -t273, 0, -pkin(7) * t272 + t227 * t256, pkin(7) * t273 + t229 * t256, -t124 * t195 + t125 * t193 - t129 * t183 + t130 * t181 - t140 * t185 - t141 * t182 - t150 * t173 - t151 * t236, t124 * t150 + t125 * t151 + t140 * t129 + t141 * t130 + (t203 + t251) * t216, -t100 * t242 + t143 * t95, t100 * t138 + t101 * t242 + t143 * t233 + t241 * t95, t100 * t218, -t101 * t218, 0, t101 * t146 - t138 * t156 - t147 * t241 - t165 * t233 - t218 * t45, t100 * t146 + t143 * t147 - t156 * t242 + t165 * t95 - t218 * t44, t101 * t38 - t11 * t138 + t120 * t45 + t221 * t237 - t233 * t53 - t241 * t8, -t101 * t39 + t12 * t138 + t122 * t45 + t223 * t237 + t233 * t54 + t241 * t9, -t11 * t122 - t12 * t120 + (-t100 * t38 - t143 * t8 - t53 * t95) * t223 + (-t100 * t39 - t143 * t9 - t54 * t95) * t221, t11 * t38 + t12 * t39 - t29 * t309 + t45 * t63 + t53 * t8 + t54 * t9, -t22 * t99 - t246 * t35, -t22 * t98 + t23 * t99 + t246 * t36 - t35 * t74, -t101 * t246 + t134 * t35 - t22 * t241 + t233 * t99, -t101 * t74 - t134 * t36 + t23 * t241 + t233 * t98, t101 * t134 + t233 * t241 (-t225 * t10 + t228 * t6) * t134 - t249 * t233 - t263 * t241 + t3 * t101 + t21 * t74 + t58 * t23 + t15 * t98 + t55 * t36 + (-t134 * t248 + t241 * t4) * qJD(6) -(t228 * t10 + t225 * t6) * t134 + t248 * t233 + t253 * t241 - t4 * t101 - t21 * t246 + t58 * t22 - t15 * t99 + t55 * t35 + (-t134 * t249 + t241 * t3) * qJD(6); 0, 0, 0, -t227 * t274, t271 * t231, 0, 0, 0, t231 * pkin(1) * t227, pkin(1) * t274 (t141 + t144) * t183 + (-t145 + t140) * t181 + (-t224 * t173 - t222 * t236) * pkin(2), -t140 * t144 - t141 * t145 + (t124 * t224 + t125 * t222 - t203 * t270) * pkin(2), t321, t320, t329, t339, 0, t138 * t155 - t218 * t289 + t315, t155 * t242 + (-t245 + t78) * t218 + t313, t120 * t289 + t138 * t42 + t221 * t234 + t314, t122 * t289 - t138 * t43 + t223 * t234 + t318, t122 * t42 + t311 * t120 + (t122 * t161 + t319) * t221 + t337, -t161 * t250 + t176 * t252 + t177 * t29 + t289 * t63 - t38 * t42 - t39 * t43, t333, t331, t332, t338, t328 -(t148 * t228 - t149 * t225) * t233 + t157 * t23 + t290 * t74 + (t225 * t244 - t228 * t243) * t134 + t334 (t148 * t225 + t149 * t228) * t233 + t157 * t22 - t290 * t246 + (t225 * t243 + t228 * t244) * t134 + t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181 ^ 2 - t183 ^ 2, t140 * t183 - t141 * t181 + t211, 0, 0, 0, 0, 0, -t233 - t286, t95 + t283, t120 * t242 - t138 * t326 - t223 * t233, t122 * t242 - t138 * t325 + t221 * t233 (-t221 ^ 2 - t223 ^ 2) * t95 + (t120 * t223 - t122 * t221) * t138, t138 * t250 + t221 * t9 + t223 * t8 + t242 * t63, 0, 0, 0, 0, 0, t254 + t297, -t293 - t344; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t321, t320, t329, t339, 0, t218 * t67 + t315, t66 * t218 + t313, -t120 * t67 + t138 * t49 + t221 * t235 + t314, -t122 * t67 - t138 * t50 + t223 * t235 + t318, t122 * t49 + t310 * t120 + (qJD(5) * t122 + t319) * t221 + t337, -pkin(4) * t29 + qJ(5) * t252 - qJD(5) * t250 - t38 * t49 - t39 * t50 - t63 * t67, t333, t331, t332, t338, t328 -(t204 * t228 - t205 * t225) * t233 + t214 * t23 - t56 * t74 + (t225 * t240 - t228 * t239) * t134 + t334 (t204 * t225 + t205 * t228) * t233 + t214 * t22 + t56 * t246 + (t225 * t239 + t228 * t240) * t134 + t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122 * t138 + t296, t120 * t138 + t294, -t120 ^ 2 - t122 ^ 2, t120 * t39 + t122 * t38 + t29, 0, 0, 0, 0, 0, t23 - t327, -t134 * t322 + (-t122 * t134 + t247) * t225 + t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t246 * t74, t246 ^ 2 - t74 ^ 2, t134 * t74 + t22, -t23 - t327, -t233, t246 * t55 + t308 * t4 + t263, t3 * t308 + t55 * t74 - t253;];
tauc_reg  = t1;
