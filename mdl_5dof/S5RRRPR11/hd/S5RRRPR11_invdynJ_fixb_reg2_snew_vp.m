% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPR11_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR11_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:15
% EndTime: 2019-12-31 21:35:26
% DurationCPUTime: 4.09s
% Computational Cost: add. (10642->358), mult. (21181->447), div. (0->0), fcn. (14029->8), ass. (0->226)
t194 = cos(qJ(3));
t190 = sin(qJ(3));
t191 = sin(qJ(2));
t233 = qJD(1) * qJD(2);
t178 = t191 * t233;
t195 = cos(qJ(2));
t232 = t195 * qJDD(1);
t166 = -t178 + t232;
t159 = -qJDD(3) + t166;
t236 = qJD(1) * t191;
t160 = -t194 * qJD(2) + t190 * t236;
t162 = t190 * qJD(2) + t194 * t236;
t251 = t162 * t160;
t201 = t159 - t251;
t245 = t190 * t201;
t158 = t162 ^ 2;
t176 = t195 * qJD(1) - qJD(3);
t272 = t176 ^ 2;
t283 = -t158 - t272;
t73 = -t194 * t283 - t245;
t318 = pkin(1) * t73;
t317 = pkin(2) * t73;
t316 = pkin(7) * t73;
t240 = t194 * t201;
t80 = -t190 * t283 + t240;
t315 = pkin(7) * t80;
t314 = t195 * t80;
t273 = t160 ^ 2;
t139 = t273 - t272;
t179 = t191 * qJDD(1);
t228 = t195 * t233;
t165 = t179 + t228;
t223 = t194 * qJDD(2) - t190 * t165;
t119 = t162 * qJD(3) - t223;
t146 = t162 * t176;
t91 = t119 + t146;
t313 = t191 * (t194 * t139 + t245) + t195 * t91;
t207 = -t190 * qJDD(2) - t194 * t165;
t120 = -t160 * qJD(3) - t207;
t252 = t160 * t176;
t284 = t120 + t252;
t261 = t190 * t284;
t282 = t158 - t273;
t286 = t119 - t146;
t312 = t191 * (t194 * t286 + t261) + t195 * t282;
t110 = t159 + t251;
t239 = t194 * t110;
t280 = -t272 - t273;
t291 = t190 * t280 - t239;
t311 = pkin(2) * t291;
t244 = t190 * t110;
t290 = t194 * t280 + t244;
t310 = pkin(7) * t290;
t309 = pkin(7) * t291;
t308 = qJ(4) * t284;
t140 = -t158 + t272;
t307 = t194 * t140 - t244;
t305 = t190 * t139 - t240;
t285 = t120 - t252;
t303 = t191 * (-t190 * t140 - t239) - t195 * t285;
t302 = pkin(6) * (t191 * t286 + t195 * t290) - pkin(1) * t291;
t281 = t158 + t273;
t301 = pkin(2) * t281;
t189 = sin(qJ(5));
t153 = qJDD(5) + t159;
t193 = cos(qJ(5));
t125 = -t193 * t160 + t189 * t162;
t127 = t189 * t160 + t193 * t162;
t83 = t127 * t125;
t287 = t153 - t83;
t300 = t189 * t287;
t298 = t191 * t281;
t296 = t193 * t287;
t197 = qJD(1) ^ 2;
t192 = sin(qJ(1));
t196 = cos(qJ(1));
t214 = t196 * g(1) + t192 * g(2);
t253 = qJDD(1) * pkin(6);
t151 = -t197 * pkin(1) - t214 + t253;
t218 = -t195 * pkin(2) - t191 * pkin(7);
t221 = t197 * t218 + t151;
t266 = t195 * g(3);
t271 = qJD(2) ^ 2;
t100 = -qJDD(2) * pkin(2) - t271 * pkin(7) + t221 * t191 + t266;
t292 = t119 * pkin(3) + t100 - t308;
t289 = -t190 * t286 + t194 * t284;
t173 = qJD(5) + t176;
t105 = t173 * t125;
t63 = -t125 * qJD(5) + t189 * t119 + t193 * t120;
t288 = -t105 + t63;
t128 = t160 * pkin(3) - t162 * qJ(4);
t267 = t191 * g(3);
t101 = -t271 * pkin(2) + qJDD(2) * pkin(7) + t221 * t195 - t267;
t226 = t192 * g(1) - t196 * g(2);
t150 = qJDD(1) * pkin(1) + t197 * pkin(6) + t226;
t210 = -t166 + t178;
t211 = t165 + t228;
t88 = t210 * pkin(2) - t211 * pkin(7) - t150;
t66 = t194 * t101 + t190 * t88;
t222 = -t159 * qJ(4) - t160 * t128 + t66;
t278 = -pkin(3) * (t283 + t272) - qJ(4) * t201 + t222;
t270 = pkin(3) + pkin(4);
t65 = t190 * t101 - t194 * t88;
t39 = t159 * pkin(3) - qJ(4) * t272 + t162 * t128 + qJDD(4) + t65;
t28 = t110 * pkin(4) - t285 * pkin(8) + t39;
t215 = t176 * pkin(4) - t162 * pkin(8);
t235 = qJD(4) * t176;
t170 = -0.2e1 * t235;
t212 = t170 + t222;
t37 = -pkin(3) * t272 + t212;
t29 = -pkin(4) * t273 + t119 * pkin(8) - t176 * t215 + t37;
t15 = t189 * t29 - t193 * t28;
t16 = t189 * t28 + t193 * t29;
t7 = -t193 * t15 + t189 * t16;
t8 = t189 * t15 + t193 * t16;
t277 = qJ(4) * t8 - t270 * t7;
t71 = t153 + t83;
t263 = t189 * t71;
t124 = t127 ^ 2;
t171 = t173 ^ 2;
t98 = -t124 - t171;
t52 = t193 * t98 - t263;
t258 = t193 * t71;
t53 = -t189 * t98 - t258;
t276 = qJ(4) * t53 - t270 * t52 + t16;
t123 = t125 ^ 2;
t77 = -t171 - t123;
t40 = t189 * t77 + t296;
t41 = t193 * t77 - t300;
t275 = qJ(4) * t41 - t270 * t40 + t15;
t225 = -t193 * t119 + t189 * t120;
t47 = (qJD(5) - t173) * t127 + t225;
t50 = t105 + t63;
t22 = -t189 * t47 - t193 * t50;
t24 = t189 * t50 - t193 * t47;
t274 = qJ(4) * t24 - t270 * t22;
t269 = pkin(3) * t194;
t38 = (-pkin(3) * t176 - 0.2e1 * qJD(4)) * t162 + t292;
t30 = t119 * pkin(4) + pkin(8) * t273 - t162 * t215 + t38;
t264 = t189 * t30;
t260 = t190 * t285;
t259 = t193 * t30;
t256 = t194 * t285;
t254 = qJ(4) * t194;
t250 = t173 * t189;
t249 = t173 * t193;
t248 = t176 * t190;
t247 = t176 * t194;
t246 = t190 * t100;
t175 = t195 * t197 * t191;
t242 = t191 * (qJDD(2) + t175);
t241 = t194 * t100;
t238 = t195 * (-t175 + qJDD(2));
t231 = t160 * t247;
t230 = t195 * t83;
t229 = t195 * t251;
t227 = -qJ(4) * t190 - pkin(2);
t35 = t190 * t65 + t194 * t66;
t136 = t191 * t151 + t266;
t137 = t195 * t151 - t267;
t224 = t191 * t136 + t195 * t137;
t138 = t162 * t248;
t219 = t191 * (t194 * t120 + t138) - t229;
t217 = -pkin(3) * t39 + qJ(4) * t37;
t216 = -pkin(3) * t285 - qJ(4) * t91;
t213 = -t194 * t119 - t160 * t248;
t209 = t190 * t66 - t194 * t65;
t208 = -pkin(1) + t218;
t204 = (t160 * t190 + t162 * t194) * t176;
t202 = t191 * (-t138 + t231) + t195 * t159;
t200 = t191 * (t190 * t119 - t231) + t229;
t199 = 0.2e1 * qJD(4) * t162 - t292;
t198 = -pkin(3) * t110 + qJ(4) * t280 - t39;
t186 = t195 ^ 2;
t185 = t191 ^ 2;
t183 = t186 * t197;
t181 = t185 * t197;
t167 = -0.2e1 * t178 + t232;
t164 = t179 + 0.2e1 * t228;
t103 = -t124 + t171;
t102 = t123 - t171;
t97 = (qJD(3) - t176) * t160 + t207;
t92 = (-qJD(3) - t176) * t162 + t223;
t85 = t190 * t120 - t162 * t247;
t82 = t124 - t123;
t69 = (-t125 * t193 + t127 * t189) * t173;
t68 = (t125 * t189 + t127 * t193) * t173;
t67 = -t123 - t124;
t62 = -t127 * qJD(5) - t225;
t61 = t194 * t92 + t260;
t60 = -t194 * t91 + t260;
t58 = -t190 * t91 - t256;
t57 = t193 * t102 - t263;
t56 = -t189 * t103 + t296;
t55 = -t189 * t102 - t258;
t54 = -t193 * t103 - t300;
t46 = (qJD(5) + t173) * t127 + t225;
t45 = -t127 * t250 + t193 * t63;
t44 = -t127 * t249 - t189 * t63;
t43 = t125 * t249 - t189 * t62;
t42 = -t125 * t250 - t193 * t62;
t36 = qJ(4) * t281 + t39;
t33 = (t281 - t272) * pkin(3) + t212;
t32 = (-t286 + t146) * pkin(3) + t199;
t31 = pkin(3) * t146 + t199 + t308;
t26 = t190 * t52 + t194 * t53;
t25 = t190 * t53 - t194 * t52;
t23 = -t189 * t288 - t193 * t46;
t21 = t189 * t46 - t193 * t288;
t20 = t190 * t40 + t194 * t41;
t19 = t190 * t41 - t194 * t40;
t18 = t190 * t39 + t194 * t37;
t17 = t190 * t37 - t194 * t39;
t14 = -pkin(8) * t52 + qJ(4) * t288 - t259;
t13 = -pkin(8) * t40 + qJ(4) * t46 - t264;
t12 = t190 * t22 + t194 * t24;
t11 = t190 * t24 - t194 * t22;
t10 = -pkin(8) * t53 + t270 * t288 + t264;
t9 = -pkin(8) * t41 + t270 * t46 - t259;
t6 = -pkin(8) * t7 - qJ(4) * t30;
t5 = -pkin(8) * t22 + qJ(4) * t67 - t7;
t4 = -pkin(8) * t24 + t270 * t67 - t8;
t3 = -pkin(8) * t8 - t270 * t30;
t2 = t190 * t7 + t194 * t8;
t1 = t190 * t8 - t194 * t7;
t27 = [0, 0, 0, 0, 0, qJDD(1), t226, t214, 0, 0, t211 * t191, t195 * t164 + t191 * t167, t242 + t195 * (-t181 + t271), -t210 * t195, t191 * (t183 - t271) + t238, 0, t195 * t150 + pkin(1) * t167 + pkin(6) * (t195 * (-t183 - t271) - t242), -t191 * t150 - pkin(1) * t164 + pkin(6) * (-t238 - t191 * (-t181 - t271)), pkin(1) * (t181 + t183) + (t185 + t186) * t253 + t224, pkin(1) * t150 + pkin(6) * t224, t219, -t312, t303, t200, t313, t202, t191 * (t246 - t309) + t195 * (t65 - t311) + t302, t191 * (t241 + t316) + t195 * (t66 + t317) + t318 + pkin(6) * (-t191 * t97 + t314), -t191 * t209 + pkin(6) * (t195 * t61 - t298) + t208 * (t190 * t92 - t256), pkin(6) * (t191 * t100 + t195 * t35) + t208 * t209, t219, t303, t312, t202, -t313, t200, t191 * (-t190 * t32 - t254 * t286 - t309) + t195 * (-t198 - t311) + t302, t191 * (-pkin(7) * t58 - t190 * t33 + t194 * t36) + t195 * (-pkin(2) * t58 - t216) - pkin(1) * t58 + pkin(6) * (t195 * t60 - t298), t191 * (-pkin(3) * t261 + t194 * t31 - t316) + t195 * (0.2e1 * t235 - t278 - t317) - t318 + pkin(6) * (-t191 * t284 - t314), t191 * (-pkin(7) * t17 + (pkin(3) * t190 - t254) * t38) + t195 * (-pkin(2) * t17 - t217) - pkin(1) * t17 + pkin(6) * (t195 * t18 + t191 * t38), t191 * (-t190 * t44 + t194 * t45) + t230, t191 * (-t190 * t21 + t194 * t23) + t195 * t82, t191 * (-t190 * t54 + t194 * t56) + t195 * t50, t191 * (-t190 * t42 + t194 * t43) - t230, t191 * (-t190 * t55 + t194 * t57) - t195 * t47, t191 * (-t190 * t68 + t194 * t69) + t195 * t153, t191 * (-pkin(7) * t19 + t194 * t13 - t190 * t9) + t195 * (-pkin(2) * t19 - t275) - pkin(1) * t19 + pkin(6) * (-t191 * t46 + t195 * t20), t191 * (-pkin(7) * t25 - t190 * t10 + t194 * t14) + t195 * (-pkin(2) * t25 - t276) - pkin(1) * t25 + pkin(6) * (-t191 * t288 + t195 * t26), t191 * (-pkin(7) * t11 - t190 * t4 + t194 * t5) + t195 * (-pkin(2) * t11 - t274) - pkin(1) * t11 + pkin(6) * (t195 * t12 - t191 * t67), t191 * (-pkin(7) * t1 - t190 * t3 + t194 * t6) + t195 * (-pkin(2) * t1 - t277) - pkin(1) * t1 + pkin(6) * (t191 * t30 + t195 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, t181 - t183, t179, t175, t232, qJDD(2), -t136, -t137, 0, 0, t85, t289, t307, t213, t305, t204, -pkin(2) * t286 - t241 + t310, pkin(2) * t97 + t246 + t315, pkin(7) * t61 + t301 + t35, -pkin(2) * t100 + pkin(7) * t35, t85, t307, -t289, t204, -t305, t213, t194 * t32 + t227 * t286 + t310, pkin(7) * t60 + t190 * t36 + t194 * t33 + t301, -t315 + t190 * t31 + (pkin(2) + t269) * t284, pkin(7) * t18 + (t227 - t269) * t38, t190 * t45 + t194 * t44, t190 * t23 + t194 * t21, t190 * t56 + t194 * t54, t190 * t43 + t194 * t42, t190 * t57 + t194 * t55, t190 * t69 + t194 * t68, pkin(2) * t46 + pkin(7) * t20 + t190 * t13 + t194 * t9, pkin(2) * t288 + pkin(7) * t26 + t194 * t10 + t190 * t14, pkin(2) * t67 + pkin(7) * t12 + t190 * t5 + t194 * t4, -pkin(2) * t30 + pkin(7) * t2 + t190 * t6 + t194 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t251, t282, t285, -t251, -t91, -t159, -t65, -t66, 0, 0, t251, t285, -t282, -t159, t91, -t251, t198, t216, t170 + t278, t217, -t83, -t82, -t50, t83, t47, -t153, t275, t276, t274, t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t285, t283, t39, 0, 0, 0, 0, 0, 0, t40, t52, t22, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t82, t50, -t83, -t47, t153, -t15, -t16, 0, 0;];
tauJ_reg = t27;
