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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:34:09
% EndTime: 2019-12-05 18:34:15
% DurationCPUTime: 3.06s
% Computational Cost: add. (5979->341), mult. (9185->418), div. (0->0), fcn. (6642->16), ass. (0->207)
t205 = cos(qJ(2));
t262 = qJD(1) * t205;
t251 = pkin(1) * t262;
t228 = qJD(3) - t251;
t197 = cos(pkin(9));
t204 = cos(qJ(4));
t267 = t204 * t197;
t196 = sin(pkin(9));
t200 = sin(qJ(4));
t270 = t200 * t196;
t137 = -t267 + t270;
t198 = -pkin(7) - qJ(3);
t151 = t198 * t196;
t185 = t197 * pkin(7);
t152 = t197 * qJ(3) + t185;
t258 = qJD(4) * t204;
t277 = t151 * t258 + qJD(3) * t267 + (-qJD(3) * t196 - qJD(4) * t152) * t200 + t137 * t251;
t138 = t204 * t196 + t200 * t197;
t96 = t200 * t151 + t204 * t152;
t276 = -qJD(4) * t96 - t228 * t138;
t199 = sin(qJ(5));
t203 = cos(qJ(5));
t194 = qJD(1) + qJD(2);
t250 = t194 * t270;
t113 = -t194 * t267 + t250;
t115 = t138 * t194;
t222 = t199 * t113 - t203 * t115;
t259 = qJD(4) * t200;
t247 = t196 * t259;
t188 = qJDD(1) + qJDD(2);
t246 = t197 * t258;
t249 = t138 * t188 + t194 * t246;
t71 = t194 * t247 - t249;
t128 = t138 * qJD(4);
t224 = t137 * t188;
t72 = t128 * t194 + t224;
t210 = qJD(5) * t222 + t199 * t71 - t203 * t72;
t193 = qJD(4) + qJD(5);
t278 = t222 * t193;
t303 = t210 - t278;
t256 = qJD(5) * t203;
t257 = qJD(5) * t199;
t217 = -t113 * t256 - t115 * t257 - t199 * t72 - t203 * t71;
t64 = -t203 * t113 - t199 * t115;
t281 = t193 * t64;
t302 = t217 - t281;
t285 = t64 ^ 2;
t286 = t222 ^ 2;
t301 = -t285 + t286;
t284 = t64 * t222;
t127 = -t246 + t247;
t292 = t127 * pkin(8);
t300 = t292 + t276;
t125 = t128 * pkin(8);
t299 = -t125 + t277;
t201 = sin(qJ(2));
t255 = qJDD(1) * t201;
t260 = qJD(2) * t205;
t100 = t188 * qJ(3) + t194 * qJD(3) + (qJD(1) * t260 + t255) * pkin(1);
t190 = t196 ^ 2;
t191 = t197 ^ 2;
t263 = t190 + t191;
t238 = t263 * t100;
t195 = qJ(1) + qJ(2);
t183 = sin(t195);
t184 = cos(t195);
t227 = g(2) * t184 + g(3) * t183;
t176 = g(2) * t183;
t265 = -g(3) * t184 + t176;
t192 = pkin(9) + qJ(4);
t182 = qJ(5) + t192;
t170 = sin(t182);
t171 = cos(t182);
t273 = t171 * t184;
t274 = t171 * t183;
t244 = pkin(7) * t188 + t100;
t81 = t244 * t196;
t82 = t244 * t197;
t240 = -t200 * t82 - t204 * t81;
t288 = t201 * pkin(1);
t252 = qJD(1) * t288;
t145 = t194 * qJ(3) + t252;
t243 = pkin(7) * t194 + t145;
t101 = t243 * t196;
t102 = t243 * t197;
t56 = -t200 * t101 + t204 * t102;
t27 = -t56 * qJD(4) + t240;
t16 = qJDD(4) * pkin(4) + t71 * pkin(8) + t27;
t254 = t101 * t258 + t200 * t81 - t204 * t82;
t26 = -t102 * t259 - t254;
t17 = -t72 * pkin(8) + t26;
t271 = t200 * t102;
t55 = -t204 * t101 - t271;
t44 = -t115 * pkin(8) + t55;
t43 = qJD(4) * pkin(4) + t44;
t45 = -t113 * pkin(8) + t56;
t4 = (qJD(5) * t43 + t17) * t203 + t199 * t16 - t45 * t257;
t173 = t197 * pkin(3) + pkin(2);
t111 = -t173 * t194 + t228;
t75 = t113 * pkin(4) + t111;
t298 = g(1) * t170 - g(2) * t274 + g(3) * t273 - t75 * t64 - t4;
t279 = t203 * t45;
t19 = t199 * t43 + t279;
t5 = -qJD(5) * t19 + t203 * t16 - t199 * t17;
t297 = -g(1) * t171 - t265 * t170 + t75 * t222 + t5;
t287 = t205 * pkin(1);
t264 = -qJD(2) * t252 + qJDD(1) * t287;
t245 = qJDD(3) - t264;
t289 = t188 * pkin(2);
t112 = t245 - t289;
t296 = -t112 + t227;
t172 = qJ(3) + t288;
t129 = (-pkin(7) - t172) * t196;
t130 = t197 * t172 + t185;
t86 = t200 * t129 + t204 * t130;
t295 = t115 ^ 2;
t202 = sin(qJ(1));
t294 = pkin(1) * t202;
t206 = cos(qJ(1));
t293 = pkin(1) * t206;
t291 = t128 * pkin(4);
t290 = t138 * pkin(8);
t95 = t204 * t151 - t200 * t152;
t76 = t95 - t290;
t132 = t137 * pkin(8);
t77 = -t132 + t96;
t38 = -t199 * t77 + t203 * t76;
t283 = qJD(5) * t38 + t300 * t199 + t299 * t203;
t39 = t199 * t76 + t203 * t77;
t282 = -qJD(5) * t39 - t299 * t199 + t300 * t203;
t280 = t199 * t45;
t275 = t115 * t113;
t272 = t197 * t188;
t266 = t227 * t197;
t261 = qJD(2) * t201;
t253 = pkin(1) * t261;
t248 = t194 * t261;
t161 = pkin(1) * t260 + qJD(3);
t237 = t263 * t161;
t236 = t263 * t188;
t85 = t204 * t129 - t200 * t130;
t181 = cos(t192);
t140 = pkin(4) * t181 + t173;
t189 = -pkin(8) + t198;
t235 = -t140 * t184 + t183 * t189;
t234 = -t173 * t184 + t183 * t198;
t90 = -t199 * t137 + t203 * t138;
t47 = qJD(5) * t90 - t199 * t127 + t203 * t128;
t91 = -t173 * t188 + t245;
t48 = t72 * pkin(4) + t91;
t89 = t203 * t137 + t199 * t138;
t233 = g(2) * t273 + g(3) * t274 + t75 * t47 + t48 * t89;
t232 = t111 * t128 + t91 * t137 + t227 * t181;
t231 = t265 + t238;
t230 = t194 * t252;
t229 = t264 + t227;
t225 = g(2) * t206 + g(3) * t202;
t59 = t85 - t290;
t60 = -t132 + t86;
t32 = -t199 * t60 + t203 * t59;
t33 = t199 * t59 + t203 * t60;
t18 = t203 * t43 - t280;
t46 = t203 * t127 + t199 * t128 + t137 * t256 + t138 * t257;
t223 = t18 * t46 - t19 * t47 - t4 * t89 - t5 * t90 + t265;
t221 = -t140 * t183 - t184 * t189;
t220 = -t173 * t183 - t184 * t198;
t109 = t137 * pkin(4) - t173;
t219 = t55 * t127 - t56 * t128 - t26 * t137 - t27 * t138 + t265;
t218 = -t252 + t291;
t216 = t230 + t289;
t178 = -pkin(2) - t287;
t215 = -pkin(1) * t248 - t178 * t188;
t214 = -t170 * t227 - t75 * t46 + t48 * t90;
t180 = sin(t192);
t213 = -t111 * t127 + t91 * t138 - t180 * t227;
t212 = -g(1) * t181 - t180 * t265;
t51 = t129 * t258 + t161 * t267 + (-qJD(4) * t130 - t161 * t196) * t200;
t211 = t228 * t263;
t52 = -qJD(4) * t86 - t138 * t161;
t187 = qJDD(4) + qJDD(5);
t168 = t184 * qJ(3);
t166 = t191 * t188;
t165 = t190 * t188;
t150 = -t173 - t287;
t146 = 0.2e1 * t196 * t272;
t139 = -t194 * pkin(2) + t228;
t110 = t113 ^ 2;
t104 = t253 + t291;
t103 = t112 * t196;
t98 = t109 - t287;
t88 = -t128 * qJD(4) - t137 * qJDD(4);
t87 = -t127 * qJD(4) + t138 * qJDD(4);
t41 = t52 + t292;
t40 = -t125 + t51;
t37 = t113 * t128 + t72 * t137;
t36 = -t115 * t127 - t71 * t138;
t29 = -t89 * t187 - t47 * t193;
t28 = t90 * t187 - t46 * t193;
t21 = t203 * t44 - t280;
t20 = -t199 * t44 - t279;
t15 = t127 * t113 - t115 * t128 + t71 * t137 - t138 * t72;
t9 = -qJD(5) * t33 - t199 * t40 + t203 * t41;
t8 = qJD(5) * t32 + t199 * t41 + t203 * t40;
t7 = -t210 * t89 - t47 * t64;
t6 = t217 * t90 + t222 * t46;
t1 = t210 * t90 - t217 * t89 + t222 * t47 - t46 * t64;
t2 = [0, 0, 0, 0, 0, qJDD(1), t225, -g(2) * t202 + g(3) * t206, 0, 0, 0, 0, 0, 0, 0, t188, (t188 * t205 - t248) * pkin(1) + t229, ((-qJDD(1) - t188) * t201 + (-qJD(1) - t194) * t260) * pkin(1) - t265, 0, (t225 + (t201 ^ 2 + t205 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t165, t146, 0, t166, 0, 0, (-t112 + t215) * t197 + t266, t103 + (-t215 - t227) * t196, t172 * t236 + t194 * t237 + t231, t112 * t178 + t139 * t253 - g(2) * (-pkin(2) * t184 - qJ(3) * t183 - t293) - g(3) * (-pkin(2) * t183 + t168 - t294) + t145 * t237 + t172 * t238, t36, t15, t87, t37, t88, 0, qJD(4) * t52 + qJDD(4) * t85 + t113 * t253 + t150 * t72 + t232, -t51 * qJD(4) - t86 * qJDD(4) + t115 * t253 - t150 * t71 + t213, -t113 * t51 - t115 * t52 + t71 * t85 - t72 * t86 + t219, t26 * t86 + t56 * t51 + t27 * t85 + t55 * t52 + t91 * t150 + t111 * t253 - g(2) * (t234 - t293) - g(3) * (t220 - t294), t6, t1, t28, t7, t29, 0, -t104 * t64 + t187 * t32 + t193 * t9 - t210 * t98 + t233, -t104 * t222 - t33 * t187 - t8 * t193 + t217 * t98 + t214, t210 * t33 - t217 * t32 + t222 * t9 + t64 * t8 + t223, t4 * t33 + t19 * t8 + t5 * t32 + t18 * t9 + t48 * t98 + t75 * t104 - g(2) * (t235 - t293) - g(3) * (t221 - t294); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, t229 + t230, (-t255 + (-qJD(2) + t194) * t262) * pkin(1) - t265, 0, 0, t165, t146, 0, t166, 0, 0, (-t112 + t216) * t197 + t266, t103 + (-t216 - t227) * t196, qJ(3) * t236 + t194 * t211 + t231, -t139 * t252 - g(3) * t168 + t296 * pkin(2) + (t238 + t176) * qJ(3) + t211 * t145, t36, t15, t87, t37, t88, 0, qJD(4) * t276 + t95 * qJDD(4) - t113 * t252 - t173 * t72 + t232, -qJD(4) * t277 - t96 * qJDD(4) - t115 * t252 + t173 * t71 + t213, -t113 * t277 - t115 * t276 + t95 * t71 - t96 * t72 + t219, -g(2) * t234 - g(3) * t220 - t111 * t252 - t91 * t173 + t26 * t96 + t27 * t95 + t276 * t55 + t277 * t56, t6, t1, t28, t7, t29, 0, -t109 * t210 + t38 * t187 + t193 * t282 - t218 * t64 + t233, t109 * t217 - t39 * t187 - t193 * t283 - t218 * t222 + t214, t210 * t39 - t217 * t38 + t222 * t282 + t283 * t64 + t223, -g(2) * t235 - g(3) * t221 + t48 * t109 + t18 * t282 + t19 * t283 + t218 * t75 + t5 * t38 + t4 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t272, t196 * t188, -t263 * t194 ^ 2, -t145 * t194 * t263 - t296, 0, 0, 0, 0, 0, 0, 0.2e1 * t115 * qJD(4) + t224, (-t113 - t250) * qJD(4) + t249, -t110 - t295, t56 * t113 + t55 * t115 - t227 + t91, 0, 0, 0, 0, 0, 0, -t210 - t278, t217 + t281, -t285 - t286, -t18 * t222 - t19 * t64 - t227 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, -t110 + t295, (t113 - t250) * qJD(4) + t249, -t275, -t224, qJDD(4), -t111 * t115 + t212 + t240, g(1) * t180 + t111 * t113 - t265 * t181 + (t55 + t271) * qJD(4) + t254, 0, 0, t284, t301, t302, -t284, t303, t187, -t20 * t193 + (t115 * t64 + t187 * t203 - t193 * t257) * pkin(4) + t297, t21 * t193 + (t115 * t222 - t187 * t199 - t193 * t256) * pkin(4) + t298, t18 * t64 - t19 * t222 - t20 * t222 - t21 * t64 + (t199 * t210 - t203 * t217 + (-t199 * t222 + t203 * t64) * qJD(5)) * pkin(4), -t18 * t20 - t19 * t21 + (-t115 * t75 + t199 * t4 + t203 * t5 + (-t18 * t199 + t19 * t203) * qJD(5) + t212) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, t301, t302, -t284, t303, t187, t19 * t193 + t297, t18 * t193 + t298, 0, 0;];
tau_reg = t2;
