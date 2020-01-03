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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:04:03
% EndTime: 2020-01-03 12:04:08
% DurationCPUTime: 3.03s
% Computational Cost: add. (5979->338), mult. (9185->418), div. (0->0), fcn. (6642->16), ass. (0->207)
t213 = cos(qJ(2));
t270 = qJD(1) * t213;
t259 = pkin(1) * t270;
t235 = qJD(3) - t259;
t205 = cos(pkin(9));
t212 = cos(qJ(4));
t277 = t212 * t205;
t204 = sin(pkin(9));
t208 = sin(qJ(4));
t280 = t208 * t204;
t139 = -t277 + t280;
t206 = -pkin(7) - qJ(3);
t155 = t206 * t204;
t191 = t205 * pkin(7);
t156 = t205 * qJ(3) + t191;
t266 = qJD(4) * t212;
t287 = t155 * t266 + qJD(3) * t277 + (-qJD(3) * t204 - qJD(4) * t156) * t208 + t139 * t259;
t140 = t212 * t204 + t208 * t205;
t96 = t208 * t155 + t212 * t156;
t286 = -qJD(4) * t96 - t235 * t140;
t209 = sin(qJ(2));
t298 = t209 * pkin(1);
t260 = qJD(1) * t298;
t297 = t213 * pkin(1);
t272 = -qJD(2) * t260 + qJDD(1) * t297;
t252 = qJDD(3) - t272;
t196 = qJDD(1) + qJDD(2);
t299 = t196 * pkin(2);
t112 = t252 - t299;
t203 = qJ(1) + qJ(2);
t189 = sin(t203);
t190 = cos(t203);
t305 = g(2) * t190 + g(3) * t189;
t313 = t112 + t305;
t207 = sin(qJ(5));
t211 = cos(qJ(5));
t202 = qJD(1) + qJD(2);
t258 = t202 * t280;
t113 = -t202 * t277 + t258;
t115 = t140 * t202;
t229 = t207 * t113 - t211 * t115;
t267 = qJD(4) * t208;
t254 = t204 * t267;
t253 = t205 * t266;
t256 = t140 * t196 + t202 * t253;
t71 = t202 * t254 - t256;
t130 = t140 * qJD(4);
t231 = t139 * t196;
t72 = t202 * t130 + t231;
t218 = t229 * qJD(5) + t207 * t71 - t211 * t72;
t201 = qJD(4) + qJD(5);
t289 = t229 * t201;
t312 = t218 - t289;
t264 = qJD(5) * t211;
t265 = qJD(5) * t207;
t225 = -t113 * t264 - t115 * t265 - t207 * t72 - t211 * t71;
t64 = -t211 * t113 - t207 * t115;
t288 = t64 * t201;
t311 = t225 - t288;
t295 = t229 ^ 2;
t296 = t64 ^ 2;
t310 = t295 - t296;
t294 = t64 * t229;
t129 = -t253 + t254;
t302 = t129 * pkin(8);
t309 = t302 + t286;
t127 = t130 * pkin(8);
t308 = -t127 + t287;
t263 = qJDD(1) * t209;
t268 = qJD(2) * t213;
t100 = t196 * qJ(3) + t202 * qJD(3) + (qJD(1) * t268 + t263) * pkin(1);
t198 = t204 ^ 2;
t199 = t205 ^ 2;
t271 = t198 + t199;
t245 = t271 * t100;
t181 = g(3) * t190;
t304 = g(2) * t189 - t181;
t200 = pkin(9) + qJ(4);
t188 = qJ(5) + t200;
t174 = sin(t188);
t175 = cos(t188);
t251 = pkin(7) * t196 + t100;
t81 = t251 * t204;
t82 = t251 * t205;
t247 = -t208 * t82 - t212 * t81;
t147 = t202 * qJ(3) + t260;
t250 = pkin(7) * t202 + t147;
t101 = t250 * t204;
t102 = t250 * t205;
t56 = -t208 * t101 + t212 * t102;
t27 = -t56 * qJD(4) + t247;
t16 = qJDD(4) * pkin(4) + t71 * pkin(8) + t27;
t262 = t101 * t266 + t208 * t81 - t212 * t82;
t26 = -t102 * t267 - t262;
t17 = -t72 * pkin(8) + t26;
t281 = t208 * t102;
t55 = -t212 * t101 - t281;
t44 = -t115 * pkin(8) + t55;
t43 = qJD(4) * pkin(4) + t44;
t45 = -t113 * pkin(8) + t56;
t4 = (qJD(5) * t43 + t17) * t211 + t207 * t16 - t45 * t265;
t177 = t205 * pkin(3) + pkin(2);
t111 = -t177 * t202 + t235;
t75 = t113 * pkin(4) + t111;
t307 = g(1) * t174 + t304 * t175 - t75 * t64 - t4;
t283 = t174 * t190;
t284 = t174 * t189;
t290 = t211 * t45;
t19 = t207 * t43 + t290;
t5 = -t19 * qJD(5) + t211 * t16 - t207 * t17;
t306 = -g(1) * t175 + g(2) * t284 - g(3) * t283 + t75 * t229 + t5;
t176 = qJ(3) + t298;
t131 = (-pkin(7) - t176) * t204;
t132 = t205 * t176 + t191;
t86 = t208 * t131 + t212 * t132;
t303 = t115 ^ 2;
t301 = t130 * pkin(4);
t300 = t140 * pkin(8);
t95 = t212 * t155 - t208 * t156;
t76 = t95 - t300;
t134 = t139 * pkin(8);
t77 = -t134 + t96;
t38 = -t207 * t77 + t211 * t76;
t293 = t38 * qJD(5) + t309 * t207 + t308 * t211;
t39 = t207 * t76 + t211 * t77;
t292 = -t39 * qJD(5) - t308 * t207 + t309 * t211;
t291 = t207 * t45;
t285 = t115 * t113;
t282 = t205 * t196;
t187 = cos(t200);
t142 = pkin(4) * t187 + t177;
t197 = -pkin(8) + t206;
t276 = t189 * t142 + t190 * t197;
t275 = t189 * t177 + t190 * t206;
t274 = t190 * pkin(2) + t189 * qJ(3);
t269 = qJD(2) * t209;
t261 = pkin(1) * t269;
t257 = t313 * t204;
t255 = t202 * t269;
t165 = pkin(1) * t268 + qJD(3);
t244 = t271 * t165;
t243 = t271 * t196;
t85 = t212 * t131 - t208 * t132;
t242 = t190 * t142 - t189 * t197;
t241 = t190 * t177 - t189 * t206;
t46 = t211 * t129 + t207 * t130 + t139 * t264 + t140 * t265;
t91 = -t177 * t196 + t252;
t48 = t72 * pkin(4) + t91;
t90 = -t207 * t139 + t211 * t140;
t240 = g(2) * t283 + g(3) * t284 - t75 * t46 + t48 * t90;
t186 = sin(t200);
t239 = -t111 * t129 + t91 * t140 + t305 * t186;
t238 = -t304 + t245;
t237 = t202 * t260;
t236 = t272 - t305;
t210 = sin(qJ(1));
t214 = cos(qJ(1));
t232 = -g(2) * t214 - g(3) * t210;
t59 = t85 - t300;
t60 = -t134 + t86;
t32 = -t207 * t60 + t211 * t59;
t33 = t207 * t59 + t211 * t60;
t18 = t211 * t43 - t291;
t47 = t90 * qJD(5) - t207 * t129 + t211 * t130;
t89 = t211 * t139 + t207 * t140;
t230 = t18 * t46 - t19 * t47 - t4 * t89 - t5 * t90 - t304;
t109 = t139 * pkin(4) - t177;
t228 = t55 * t129 - t56 * t130 - t26 * t139 - t27 * t140 - t304;
t226 = -t260 + t301;
t224 = -t237 - t299;
t184 = -pkin(2) - t297;
t223 = pkin(1) * t255 + t184 * t196;
t222 = -t175 * t305 + t75 * t47 + t48 * t89;
t221 = t111 * t130 + t91 * t139 - t187 * t305;
t220 = -g(1) * t187 + t186 * t304;
t51 = t131 * t266 + t165 * t277 + (-qJD(4) * t132 - t165 * t204) * t208;
t219 = t235 * t271;
t52 = -qJD(4) * t86 - t140 * t165;
t195 = qJDD(4) + qJDD(5);
t193 = t214 * pkin(1);
t192 = t210 * pkin(1);
t178 = t189 * pkin(2);
t170 = t199 * t196;
t169 = t198 * t196;
t154 = -t177 - t297;
t150 = 0.2e1 * t204 * t282;
t141 = -t202 * pkin(2) + t235;
t110 = t113 ^ 2;
t104 = t261 + t301;
t98 = t109 - t297;
t88 = -t130 * qJD(4) - t139 * qJDD(4);
t87 = -t129 * qJD(4) + t140 * qJDD(4);
t41 = t52 + t302;
t40 = -t127 + t51;
t37 = t113 * t130 + t72 * t139;
t36 = -t115 * t129 - t71 * t140;
t29 = -t89 * t195 - t47 * t201;
t28 = t90 * t195 - t46 * t201;
t21 = t211 * t44 - t291;
t20 = -t207 * t44 - t290;
t15 = t129 * t113 - t115 * t130 + t71 * t139 - t140 * t72;
t9 = -t33 * qJD(5) - t207 * t40 + t211 * t41;
t8 = t32 * qJD(5) + t207 * t41 + t211 * t40;
t7 = -t218 * t89 - t47 * t64;
t6 = t225 * t90 + t229 * t46;
t1 = t218 * t90 - t225 * t89 + t229 * t47 - t46 * t64;
t2 = [0, 0, 0, 0, 0, qJDD(1), t232, g(2) * t210 - g(3) * t214, 0, 0, 0, 0, 0, 0, 0, t196, (t196 * t213 - t255) * pkin(1) + t236, ((-qJDD(1) - t196) * t209 + (-qJD(1) - t202) * t268) * pkin(1) + t304, 0, (t232 + (t209 ^ 2 + t213 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t169, t150, 0, t170, 0, 0, (-t223 - t313) * t205, t204 * t223 + t257, t176 * t243 + t202 * t244 + t238, t112 * t184 + t141 * t261 - g(2) * (t193 + t274) - g(3) * (-t190 * qJ(3) + t178 + t192) + t147 * t244 + t176 * t245, t36, t15, t87, t37, t88, 0, t52 * qJD(4) + t85 * qJDD(4) + t113 * t261 + t154 * t72 + t221, -t51 * qJD(4) - t86 * qJDD(4) + t115 * t261 - t154 * t71 + t239, -t51 * t113 - t52 * t115 + t85 * t71 - t86 * t72 + t228, t26 * t86 + t56 * t51 + t27 * t85 + t55 * t52 + t91 * t154 + t111 * t261 - g(2) * (t193 + t241) - g(3) * (t192 + t275), t6, t1, t28, t7, t29, 0, -t104 * t64 + t32 * t195 + t9 * t201 - t218 * t98 + t222, -t104 * t229 - t33 * t195 - t8 * t201 + t225 * t98 + t240, t218 * t33 - t225 * t32 + t229 * t9 + t64 * t8 + t230, t4 * t33 + t19 * t8 + t5 * t32 + t18 * t9 + t48 * t98 + t75 * t104 - g(2) * (t193 + t242) - g(3) * (t192 + t276); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, t236 + t237, (-t263 + (-qJD(2) + t202) * t270) * pkin(1) + t304, 0, 0, t169, t150, 0, t170, 0, 0, (-t224 - t313) * t205, t204 * t224 + t257, qJ(3) * t243 + t202 * t219 + t238, -t112 * pkin(2) - t141 * t260 - g(2) * t274 - g(3) * t178 + (t245 + t181) * qJ(3) + t219 * t147, t36, t15, t87, t37, t88, 0, qJD(4) * t286 + t95 * qJDD(4) - t113 * t260 - t177 * t72 + t221, -qJD(4) * t287 - t96 * qJDD(4) - t115 * t260 + t177 * t71 + t239, -t113 * t287 - t115 * t286 + t95 * t71 - t96 * t72 + t228, -g(2) * t241 - g(3) * t275 - t111 * t260 - t91 * t177 + t26 * t96 + t27 * t95 + t286 * t55 + t287 * t56, t6, t1, t28, t7, t29, 0, -t109 * t218 + t38 * t195 + t201 * t292 - t226 * t64 + t222, t109 * t225 - t39 * t195 - t201 * t293 - t226 * t229 + t240, t218 * t39 - t225 * t38 + t229 * t292 + t293 * t64 + t230, -g(2) * t242 - g(3) * t276 + t48 * t109 + t18 * t292 + t19 * t293 + t226 * t75 + t5 * t38 + t4 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t282, t204 * t196, -t271 * t202 ^ 2, -t147 * t202 * t271 + t313, 0, 0, 0, 0, 0, 0, 0.2e1 * t115 * qJD(4) + t231, (-t113 - t258) * qJD(4) + t256, -t110 - t303, t56 * t113 + t55 * t115 + t305 + t91, 0, 0, 0, 0, 0, 0, -t218 - t289, t225 + t288, -t295 - t296, -t18 * t229 - t19 * t64 + t305 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, -t110 + t303, (t113 - t258) * qJD(4) + t256, -t285, -t231, qJDD(4), -t111 * t115 + t220 + t247, g(1) * t186 + t111 * t113 + t304 * t187 + (t55 + t281) * qJD(4) + t262, 0, 0, t294, t310, t311, -t294, t312, t195, -t20 * t201 + (t115 * t64 + t195 * t211 - t201 * t265) * pkin(4) + t306, t21 * t201 + (t115 * t229 - t195 * t207 - t201 * t264) * pkin(4) + t307, t18 * t64 - t19 * t229 - t20 * t229 - t21 * t64 + (t207 * t218 - t211 * t225 + (-t207 * t229 + t211 * t64) * qJD(5)) * pkin(4), -t18 * t20 - t19 * t21 + (-t115 * t75 + t207 * t4 + t211 * t5 + (-t18 * t207 + t19 * t211) * qJD(5) + t220) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t294, t310, t311, -t294, t312, t195, t19 * t201 + t306, t18 * t201 + t307, 0, 0;];
tau_reg = t2;
