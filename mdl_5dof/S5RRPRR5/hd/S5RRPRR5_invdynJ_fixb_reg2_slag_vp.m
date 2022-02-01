% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:55
% EndTime: 2022-01-20 11:03:07
% DurationCPUTime: 3.46s
% Computational Cost: add. (5979->344), mult. (9185->423), div. (0->0), fcn. (6642->16), ass. (0->211)
t211 = cos(qJ(2));
t269 = qJD(1) * t211;
t258 = pkin(1) * t269;
t234 = qJD(3) - t258;
t203 = cos(pkin(9));
t210 = cos(qJ(4));
t275 = t210 * t203;
t202 = sin(pkin(9));
t206 = sin(qJ(4));
t278 = t206 * t202;
t138 = -t275 + t278;
t204 = -pkin(7) - qJ(3);
t155 = t204 * t202;
t190 = t203 * pkin(7);
t156 = t203 * qJ(3) + t190;
t265 = qJD(4) * t210;
t287 = t155 * t265 + qJD(3) * t275 + (-qJD(3) * t202 - qJD(4) * t156) * t206 + t138 * t258;
t139 = t210 * t202 + t206 * t203;
t96 = t206 * t155 + t210 * t156;
t286 = -qJD(4) * t96 - t139 * t234;
t205 = sin(qJ(5));
t209 = cos(qJ(5));
t200 = qJD(1) + qJD(2);
t257 = t200 * t278;
t113 = -t200 * t275 + t257;
t115 = t139 * t200;
t229 = t205 * t113 - t209 * t115;
t266 = qJD(4) * t206;
t254 = t202 * t266;
t194 = qJDD(1) + qJDD(2);
t253 = t203 * t265;
t256 = t139 * t194 + t200 * t253;
t71 = t200 * t254 - t256;
t129 = t139 * qJD(4);
t231 = t138 * t194;
t72 = t129 * t200 + t231;
t216 = qJD(5) * t229 + t205 * t71 - t209 * t72;
t199 = qJD(4) + qJD(5);
t289 = t229 * t199;
t314 = t216 - t289;
t263 = qJD(5) * t209;
t264 = qJD(5) * t205;
t221 = -t113 * t263 - t115 * t264 - t205 * t72 - t209 * t71;
t64 = -t209 * t113 - t205 * t115;
t288 = t64 * t199;
t313 = t221 - t288;
t295 = t229 ^ 2;
t296 = t64 ^ 2;
t312 = t295 - t296;
t294 = t64 * t229;
t128 = -t253 + t254;
t303 = t128 * pkin(8);
t311 = t303 + t286;
t126 = t129 * pkin(8);
t310 = -t126 + t287;
t207 = sin(qJ(2));
t262 = qJDD(1) * t207;
t267 = qJD(2) * t211;
t100 = t194 * qJ(3) + t200 * qJD(3) + (qJD(1) * t267 + t262) * pkin(1);
t196 = t202 ^ 2;
t197 = t203 ^ 2;
t270 = t196 + t197;
t243 = t270 * t100;
t299 = t207 * pkin(1);
t259 = qJD(1) * t299;
t297 = t211 * pkin(1);
t271 = -qJD(2) * t259 + qJDD(1) * t297;
t252 = qJDD(3) - t271;
t300 = t194 * pkin(2);
t112 = t252 - t300;
t201 = qJ(1) + qJ(2);
t189 = cos(t201);
t180 = g(2) * t189;
t309 = t112 + t180;
t188 = sin(t201);
t181 = g(1) * t188;
t305 = t180 - t181;
t198 = pkin(9) + qJ(4);
t187 = qJ(5) + t198;
t174 = sin(t187);
t175 = cos(t187);
t281 = t175 * t189;
t282 = t175 * t188;
t250 = pkin(7) * t194 + t100;
t81 = t250 * t202;
t82 = t250 * t203;
t245 = -t206 * t82 - t210 * t81;
t148 = t200 * qJ(3) + t259;
t249 = pkin(7) * t200 + t148;
t101 = t249 * t202;
t102 = t249 * t203;
t56 = -t206 * t101 + t210 * t102;
t27 = -t56 * qJD(4) + t245;
t16 = qJDD(4) * pkin(4) + t71 * pkin(8) + t27;
t261 = t101 * t265 + t206 * t81 - t210 * t82;
t26 = -t102 * t266 - t261;
t17 = -t72 * pkin(8) + t26;
t279 = t206 * t102;
t55 = -t210 * t101 - t279;
t44 = -t115 * pkin(8) + t55;
t43 = qJD(4) * pkin(4) + t44;
t45 = -t113 * pkin(8) + t56;
t4 = (qJD(5) * t43 + t17) * t209 + t205 * t16 - t45 * t264;
t177 = t203 * pkin(3) + pkin(2);
t111 = -t177 * t200 + t234;
t75 = t113 * pkin(4) + t111;
t308 = g(1) * t281 + g(2) * t282 + g(3) * t174 - t75 * t64 - t4;
t283 = t174 * t189;
t284 = t174 * t188;
t290 = t209 * t45;
t19 = t205 * t43 + t290;
t5 = -qJD(5) * t19 + t209 * t16 - t205 * t17;
t307 = g(1) * t283 + g(2) * t284 - g(3) * t175 + t75 * t229 + t5;
t176 = qJ(3) + t299;
t130 = (-pkin(7) - t176) * t202;
t131 = t203 * t176 + t190;
t86 = t206 * t130 + t210 * t131;
t306 = g(1) * t189 + g(2) * t188;
t304 = t115 ^ 2;
t302 = t129 * pkin(4);
t301 = t139 * pkin(8);
t208 = sin(qJ(1));
t298 = t208 * pkin(1);
t95 = t210 * t155 - t206 * t156;
t76 = t95 - t301;
t133 = t138 * pkin(8);
t77 = -t133 + t96;
t38 = -t205 * t77 + t209 * t76;
t293 = qJD(5) * t38 + t205 * t311 + t209 * t310;
t39 = t205 * t76 + t209 * t77;
t292 = -qJD(5) * t39 - t205 * t310 + t209 * t311;
t291 = t205 * t45;
t285 = t115 * t113;
t280 = t203 * t194;
t274 = t309 * t202;
t273 = t189 * pkin(2) + t188 * qJ(3);
t268 = qJD(2) * t207;
t260 = pkin(1) * t268;
t255 = t200 * t268;
t251 = -pkin(2) * t188 + t189 * qJ(3);
t164 = pkin(1) * t267 + qJD(3);
t242 = t270 * t164;
t241 = t270 * t194;
t85 = t210 * t130 - t206 * t131;
t186 = cos(t198);
t141 = pkin(4) * t186 + t177;
t195 = pkin(8) - t204;
t240 = t189 * t141 + t188 * t195;
t239 = -t141 * t188 + t195 * t189;
t238 = t189 * t177 - t188 * t204;
t237 = -t306 + t243;
t236 = t200 * t259;
t235 = -t271 + t305;
t212 = cos(qJ(1));
t232 = g(1) * t208 - g(2) * t212;
t59 = t85 - t301;
t60 = -t133 + t86;
t32 = -t205 * t60 + t209 * t59;
t33 = t205 * t59 + t209 * t60;
t18 = t209 * t43 - t291;
t46 = t209 * t128 + t205 * t129 + t138 * t263 + t139 * t264;
t90 = -t205 * t138 + t209 * t139;
t47 = qJD(5) * t90 - t205 * t128 + t209 * t129;
t89 = t209 * t138 + t205 * t139;
t230 = t18 * t46 - t19 * t47 - t4 * t89 - t5 * t90 - t306;
t228 = -t177 * t188 - t189 * t204;
t109 = t138 * pkin(4) - t177;
t227 = t55 * t128 - t56 * t129 - t26 * t138 - t27 * t139 - t306;
t226 = -t259 + t302;
t91 = -t177 * t194 + t252;
t48 = t72 * pkin(4) + t91;
t225 = -g(1) * t284 + g(2) * t283 - t75 * t46 + t48 * t90;
t185 = sin(t198);
t224 = -t111 * t128 + t91 * t139 + t185 * t305;
t223 = g(1) * t282 - g(2) * t281 + t75 * t47 + t48 * t89;
t222 = t111 * t129 + t91 * t138 - t186 * t305;
t220 = -t236 - t300;
t183 = -pkin(2) - t297;
t219 = pkin(1) * t255 + t183 * t194;
t218 = -g(3) * t186 + t185 * t306;
t51 = t130 * t265 + t164 * t275 + (-qJD(4) * t131 - t164 * t202) * t206;
t217 = t234 * t270;
t52 = -qJD(4) * t86 - t139 * t164;
t193 = qJDD(4) + qJDD(5);
t191 = t212 * pkin(1);
t169 = t197 * t194;
t168 = t196 * t194;
t162 = t203 * t181;
t154 = -t177 - t297;
t150 = 0.2e1 * t202 * t280;
t140 = -t200 * pkin(2) + t234;
t110 = t113 ^ 2;
t104 = t260 + t302;
t98 = t109 - t297;
t88 = -t129 * qJD(4) - t138 * qJDD(4);
t87 = -t128 * qJD(4) + t139 * qJDD(4);
t41 = t52 + t303;
t40 = -t126 + t51;
t37 = t113 * t129 + t72 * t138;
t36 = -t115 * t128 - t71 * t139;
t29 = -t89 * t193 - t47 * t199;
t28 = t90 * t193 - t46 * t199;
t21 = t209 * t44 - t291;
t20 = -t205 * t44 - t290;
t15 = t128 * t113 - t115 * t129 + t71 * t138 - t139 * t72;
t9 = -qJD(5) * t33 - t205 * t40 + t209 * t41;
t8 = qJD(5) * t32 + t205 * t41 + t209 * t40;
t7 = -t216 * t89 - t47 * t64;
t6 = t221 * t90 + t229 * t46;
t1 = t216 * t90 - t221 * t89 + t229 * t47 - t46 * t64;
t2 = [0, 0, 0, 0, 0, qJDD(1), t232, g(1) * t212 + g(2) * t208, 0, 0, 0, 0, 0, 0, 0, t194, (t194 * t211 - t255) * pkin(1) - t235, ((-qJDD(1) - t194) * t207 + (-qJD(1) - t200) * t267) * pkin(1) + t306, 0, (t232 + (t207 ^ 2 + t211 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t168, t150, 0, t169, 0, 0, t162 + (-t219 - t309) * t203, (t219 - t181) * t202 + t274, t176 * t241 + t200 * t242 + t237, t112 * t183 + t140 * t260 - g(1) * (t251 - t298) - g(2) * (t191 + t273) + t148 * t242 + t176 * t243, t36, t15, t87, t37, t88, 0, qJD(4) * t52 + qJDD(4) * t85 + t113 * t260 + t154 * t72 + t222, -qJD(4) * t51 - qJDD(4) * t86 + t115 * t260 - t154 * t71 + t224, -t113 * t51 - t115 * t52 + t71 * t85 - t72 * t86 + t227, t26 * t86 + t56 * t51 + t27 * t85 + t55 * t52 + t91 * t154 + t111 * t260 - g(1) * (t228 - t298) - g(2) * (t191 + t238), t6, t1, t28, t7, t29, 0, -t104 * t64 + t193 * t32 + t199 * t9 - t216 * t98 + t223, -t104 * t229 - t193 * t33 - t199 * t8 + t221 * t98 + t225, t216 * t33 - t221 * t32 + t229 * t9 + t64 * t8 + t230, t4 * t33 + t19 * t8 + t5 * t32 + t18 * t9 + t48 * t98 + t75 * t104 - g(1) * (t239 - t298) - g(2) * (t191 + t240); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, -t235 + t236, (-t262 + (-qJD(2) + t200) * t269) * pkin(1) + t306, 0, 0, t168, t150, 0, t169, 0, 0, t162 + (-t220 - t309) * t203, (t220 - t181) * t202 + t274, qJ(3) * t241 + t200 * t217 + t237, -t112 * pkin(2) - g(1) * t251 - g(2) * t273 + qJ(3) * t243 - t140 * t259 + t148 * t217, t36, t15, t87, t37, t88, 0, qJD(4) * t286 + t95 * qJDD(4) - t113 * t259 - t177 * t72 + t222, -qJD(4) * t287 - t96 * qJDD(4) - t115 * t259 + t177 * t71 + t224, -t113 * t287 - t115 * t286 + t95 * t71 - t96 * t72 + t227, -g(1) * t228 - g(2) * t238 - t111 * t259 - t91 * t177 + t26 * t96 + t27 * t95 + t286 * t55 + t287 * t56, t6, t1, t28, t7, t29, 0, -t109 * t216 + t38 * t193 + t199 * t292 - t226 * t64 + t223, t109 * t221 - t39 * t193 - t199 * t293 - t226 * t229 + t225, t216 * t39 - t221 * t38 + t229 * t292 + t293 * t64 + t230, -g(1) * t239 - g(2) * t240 + t48 * t109 + t18 * t292 + t19 * t293 + t226 * t75 + t5 * t38 + t4 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t280, t202 * t194, -t270 * t200 ^ 2, -t148 * t200 * t270 + t112 + t305, 0, 0, 0, 0, 0, 0, 0.2e1 * t115 * qJD(4) + t231, (-t113 - t257) * qJD(4) + t256, -t110 - t304, t56 * t113 + t55 * t115 + t305 + t91, 0, 0, 0, 0, 0, 0, -t216 - t289, t221 + t288, -t295 - t296, -t18 * t229 - t19 * t64 + t305 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, -t110 + t304, (t113 - t257) * qJD(4) + t256, -t285, -t231, qJDD(4), -t111 * t115 + t218 + t245, g(3) * t185 + t111 * t113 + t306 * t186 + (t55 + t279) * qJD(4) + t261, 0, 0, t294, t312, t313, -t294, t314, t193, -t20 * t199 + (t115 * t64 + t193 * t209 - t199 * t264) * pkin(4) + t307, t21 * t199 + (t115 * t229 - t193 * t205 - t199 * t263) * pkin(4) + t308, t18 * t64 - t19 * t229 - t20 * t229 - t21 * t64 + (t205 * t216 - t209 * t221 + (-t205 * t229 + t209 * t64) * qJD(5)) * pkin(4), -t18 * t20 - t19 * t21 + (-t115 * t75 + t205 * t4 + t209 * t5 + (-t18 * t205 + t19 * t209) * qJD(5) + t218) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t294, t312, t313, -t294, t314, t193, t19 * t199 + t307, t18 * t199 + t308, 0, 0;];
tau_reg = t2;
