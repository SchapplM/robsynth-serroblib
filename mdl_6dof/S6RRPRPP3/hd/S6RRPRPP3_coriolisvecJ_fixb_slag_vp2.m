% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:58:01
% EndTime: 2018-11-23 16:58:14
% DurationCPUTime: 13.53s
% Computational Cost: add. (6075->629), mult. (15646->809), div. (0->0), fcn. (10649->6), ass. (0->267)
t392 = Ifges(7,5) + Ifges(5,5);
t368 = Ifges(6,4) - t392;
t393 = Ifges(6,5) + Ifges(7,4);
t367 = -Ifges(5,6) + t393;
t239 = sin(qJ(4));
t241 = cos(qJ(2));
t237 = cos(pkin(9));
t332 = cos(qJ(4));
t288 = t332 * t237;
t276 = t241 * t288;
t261 = qJD(1) * t276;
t236 = sin(pkin(9));
t300 = qJD(1) * t241;
t287 = t236 * t300;
t161 = -t239 * t287 + t261;
t280 = qJD(4) * t332;
t297 = qJD(4) * t239;
t369 = -t236 * t297 + t237 * t280;
t303 = t161 - t369;
t240 = sin(qJ(2));
t301 = qJD(1) * t240;
t283 = t236 * t301;
t254 = -qJD(2) * t237 + t283;
t175 = t332 * t254;
t286 = t237 * t301;
t296 = t236 * qJD(2);
t190 = t286 + t296;
t295 = qJD(1) * qJD(2);
t277 = t241 * t295;
t274 = t236 * t277;
t85 = qJD(4) * t175 - qJD(2) * t261 + (qJD(4) * t190 + t274) * t239;
t360 = -t85 / 0.2e1;
t359 = t85 / 0.2e1;
t244 = t190 * t332 - t239 * t254;
t193 = t236 * t332 + t239 * t237;
t248 = t241 * t193;
t247 = qJD(2) * t248;
t86 = qJD(1) * t247 + qJD(4) * t244;
t357 = t86 / 0.2e1;
t278 = t240 * t295;
t397 = t278 / 0.2e1;
t385 = Ifges(5,1) + Ifges(7,3);
t384 = -Ifges(5,4) + Ifges(7,6);
t381 = Ifges(7,2) + Ifges(6,3);
t380 = Ifges(6,6) - Ifges(7,6);
t160 = qJD(1) * t248;
t178 = t193 * qJD(4);
t304 = t160 - t178;
t396 = Ifges(5,3) + Ifges(7,1) + Ifges(6,1);
t395 = t384 * t357 + t385 * t360 + t392 * t397;
t394 = t357 * t381 + t359 * t380 + t393 * t397;
t329 = pkin(8) + qJ(3);
t204 = t329 * t236;
t205 = t329 * t237;
t298 = qJD(3) * t236;
t100 = t239 * (qJD(4) * t205 + t298) - qJD(3) * t288 + t204 * t280;
t262 = pkin(2) * t240 - qJ(3) * t241;
t195 = t262 * qJD(1);
t152 = pkin(7) * t283 + t237 * t195;
t308 = t237 * t241;
t258 = pkin(3) * t240 - pkin(8) * t308;
t118 = qJD(1) * t258 + t152;
t179 = t236 * t195;
t309 = t237 * t240;
t310 = t236 * t241;
t253 = -pkin(7) * t309 - pkin(8) * t310;
t136 = qJD(1) * t253 + t179;
t55 = t239 * t118 + t332 * t136;
t41 = -qJ(5) * t301 - t55;
t391 = t100 - t41;
t150 = -t239 * t204 + t205 * t332;
t101 = qJD(3) * t193 + qJD(4) * t150;
t266 = -t118 * t332 + t239 * t136;
t390 = t101 - t266;
t232 = pkin(7) * t300;
t184 = pkin(3) * t287 + t232;
t389 = qJ(5) * t303 - qJD(5) * t193 - t184;
t200 = -pkin(2) * t241 - t240 * qJ(3) - pkin(1);
t183 = t200 * qJD(1);
t210 = qJD(2) * qJ(3) + t232;
t139 = t237 * t183 - t236 * t210;
t246 = -pkin(3) * t300 - t190 * pkin(8) + t139;
t90 = t332 * t246;
t140 = t236 * t183 + t237 * t210;
t98 = -pkin(8) * t254 + t140;
t36 = t239 * t98 - t90;
t256 = pkin(5) * t244 + t36;
t388 = qJD(5) + t256;
t133 = t190 * t239 + t175;
t387 = t133 * pkin(5) - qJD(6);
t363 = -Ifges(3,6) / 0.2e1;
t358 = -t86 / 0.2e1;
t346 = -t133 / 0.2e1;
t345 = t133 / 0.2e1;
t227 = -qJD(4) + t300;
t338 = -t227 / 0.2e1;
t344 = -t244 / 0.2e1;
t343 = t244 / 0.2e1;
t386 = -t278 / 0.2e1;
t192 = t236 * t239 - t288;
t330 = pkin(4) + qJ(6);
t377 = qJD(6) * t192 - t304 * t330 + t389;
t129 = Ifges(7,6) * t244;
t314 = t244 * Ifges(6,6);
t376 = t133 * t381 - t227 * t393 + t129 - t314;
t131 = Ifges(5,4) * t133;
t316 = Ifges(7,6) * t133;
t375 = -t227 * t392 + t244 * t385 - t131 + t316;
t374 = pkin(5) * t304 - t391;
t279 = t240 * t330;
t373 = -pkin(5) * t303 + qJD(1) * t279 + t390;
t372 = -pkin(4) * t304 + t389;
t371 = -qJD(5) - t36;
t245 = t239 * t246;
t37 = t332 * t98 + t245;
t370 = -t37 + t387;
t366 = t396 * t278 + t367 * t86 + t368 * t85;
t174 = qJD(2) * t262 - t240 * qJD(3);
t162 = t174 * qJD(1);
t231 = pkin(7) * t301;
t198 = (qJD(3) - t231) * qJD(2);
t116 = t237 * t162 - t236 * t198;
t250 = t258 * qJD(2);
t89 = qJD(1) * t250 + t116;
t117 = t236 * t162 + t237 * t198;
t97 = -pkin(8) * t274 + t117;
t255 = qJD(4) * t245 + t239 * t97 + t98 * t280 - t332 * t89;
t269 = qJD(2) * t279;
t1 = -t85 * pkin(5) - qJD(1) * t269 + t227 * qJD(6) + t255;
t6 = qJD(4) * t90 + t239 * t89 - t297 * t98 + t332 * t97;
t4 = -qJ(5) * t278 + qJD(5) * t227 - t6;
t2 = -pkin(5) * t86 - t4;
t5 = -pkin(4) * t278 + t255;
t365 = mrSges(5,1) * t255 + t6 * mrSges(5,2) - t5 * mrSges(6,2) - t2 * mrSges(7,2) + t4 * mrSges(6,3) + t1 * mrSges(7,3);
t14 = t227 * t330 + t388;
t32 = t227 * qJ(5) - t37;
t17 = -t32 - t387;
t271 = Ifges(5,6) / 0.2e1 - Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1;
t272 = Ifges(6,4) / 0.2e1 - Ifges(7,5) / 0.2e1 - Ifges(5,5) / 0.2e1;
t273 = Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1;
t31 = pkin(4) * t227 - t371;
t364 = -t271 * t133 - t273 * t227 - t272 * t244 + t139 * mrSges(4,1) + t17 * mrSges(7,2) + t31 * mrSges(6,2) - Ifges(4,6) * t254 / 0.2e1 - Ifges(4,3) * t300 / 0.2e1 + Ifges(4,5) * t190 + qJD(2) * t363 - (Ifges(3,4) * t240 + Ifges(3,2) * t241) * qJD(1) / 0.2e1 + Ifges(5,6) * t346 + Ifges(6,4) * t344 - pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t300) - t14 * mrSges(7,3) - t140 * mrSges(4,2) - t32 * mrSges(6,3) - t36 * mrSges(5,1) - t37 * mrSges(5,2) + t393 * t345 + t392 * t343 + t396 * t338;
t362 = Ifges(6,4) * t386 + Ifges(6,2) * t360 + Ifges(6,6) * t358;
t361 = Ifges(5,4) * t359 + Ifges(5,2) * t357 + Ifges(5,6) * t386;
t356 = pkin(1) * mrSges(3,1);
t355 = pkin(1) * mrSges(3,2);
t337 = t227 / 0.2e1;
t336 = -t236 / 0.2e1;
t335 = t236 / 0.2e1;
t334 = t237 / 0.2e1;
t328 = mrSges(4,2) * t237;
t327 = mrSges(5,3) * t133;
t326 = mrSges(5,3) * t244;
t325 = Ifges(4,1) * t190;
t324 = Ifges(4,1) * t237;
t323 = Ifges(4,4) * t190;
t322 = Ifges(4,4) * t236;
t321 = Ifges(4,4) * t237;
t319 = Ifges(4,5) * t237;
t318 = Ifges(4,2) * t236;
t317 = Ifges(4,6) * t236;
t315 = t244 * Ifges(5,4);
t313 = Ifges(3,5) * qJD(2);
t312 = qJ(5) * t133;
t311 = t236 * t240;
t102 = mrSges(5,2) * t227 - t327;
t105 = mrSges(6,1) * t133 + mrSges(6,3) * t227;
t307 = -t102 + t105;
t106 = -mrSges(7,1) * t133 - mrSges(7,2) * t227;
t306 = -t105 + t106;
t103 = -mrSges(5,1) * t227 - t326;
t107 = mrSges(6,1) * t244 - mrSges(6,2) * t227;
t305 = t107 - t103;
t189 = t237 * t200;
t138 = -pkin(8) * t309 + t189 + (-pkin(7) * t236 - pkin(3)) * t241;
t157 = pkin(7) * t308 + t236 * t200;
t146 = -pkin(8) * t311 + t157;
t70 = t239 * t138 + t332 * t146;
t163 = mrSges(4,1) * t274 + t277 * t328;
t299 = qJD(2) * t240;
t291 = pkin(7) * t299;
t144 = t237 * t174 + t236 * t291;
t226 = pkin(7) * t277;
t173 = pkin(3) * t274 + t226;
t285 = t241 * t296;
t185 = qJD(2) * t241 * pkin(7) + pkin(3) * t285;
t196 = pkin(3) * t311 + t240 * pkin(7);
t290 = Ifges(4,5) * t300;
t289 = Ifges(4,6) * t300;
t229 = -pkin(3) * t237 - pkin(2);
t29 = t85 * mrSges(7,2) + t86 * mrSges(7,3);
t66 = -t85 * mrSges(6,1) + mrSges(6,2) * t278;
t65 = -t86 * mrSges(7,1) + mrSges(7,2) * t278;
t149 = t332 * t204 + t205 * t239;
t199 = -qJD(2) * pkin(2) + qJD(3) + t231;
t270 = m(4) * t199 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t254 + t190 * mrSges(4,2) + mrSges(3,3) * t301;
t56 = qJ(5) * t241 - t70;
t169 = -t239 * t311 + t240 * t288;
t267 = -qJ(5) * t169 + t196;
t69 = t138 * t332 - t239 * t146;
t265 = mrSges(4,1) * t236 + t328;
t264 = -t322 + t324;
t263 = -t318 + t321;
t63 = -t85 * mrSges(7,1) - mrSges(7,3) * t278;
t257 = -qJ(5) * t193 + t229;
t57 = t241 * pkin(4) - t69;
t108 = t250 + t144;
t164 = t236 * t174;
t119 = qJD(2) * t253 + t164;
t15 = t239 * t108 + t332 * t119 + t138 * t280 - t146 * t297;
t252 = -t108 * t332 + t239 * t119 + t138 * t297 + t146 * t280;
t251 = qJ(5) * t85 - qJD(5) * t244 + t173;
t112 = -qJD(2) * t276 + t178 * t240 + t239 * t285;
t249 = qJ(5) * t112 - qJD(5) * t169 + t185;
t151 = pkin(3) * t254 + t199;
t12 = -qJ(5) * t299 + qJD(5) * t241 - t15;
t243 = -qJ(5) * t244 + t151;
t230 = Ifges(3,4) * t300;
t182 = Ifges(3,1) * t301 + t230 + t313;
t172 = (mrSges(4,1) * t240 - mrSges(4,3) * t308) * t295;
t171 = (-mrSges(4,2) * t240 - mrSges(4,3) * t310) * t295;
t168 = t193 * t240;
t159 = -mrSges(4,1) * t300 - t190 * mrSges(4,3);
t158 = mrSges(4,2) * t300 - mrSges(4,3) * t254;
t156 = -pkin(7) * t310 + t189;
t153 = -pkin(7) * t286 + t179;
t148 = (Ifges(4,5) * t240 + t241 * t264) * t295;
t147 = (Ifges(4,6) * t240 + t241 * t263) * t295;
t145 = -t237 * t291 + t164;
t130 = Ifges(6,6) * t133;
t128 = pkin(4) * t192 + t257;
t126 = -Ifges(4,4) * t254 - t290 + t325;
t125 = -Ifges(4,2) * t254 - t289 + t323;
t113 = t240 * t369 + t247;
t111 = -t192 * pkin(5) + t150;
t110 = pkin(5) * t193 + t149;
t104 = mrSges(7,1) * t244 + mrSges(7,3) * t227;
t96 = t192 * t330 + t257;
t93 = pkin(4) * t168 + t267;
t74 = t85 * mrSges(5,2);
t73 = t85 * mrSges(6,3);
t68 = -mrSges(5,2) * t278 - mrSges(5,3) * t86;
t67 = mrSges(5,1) * t278 + mrSges(5,3) * t85;
t64 = mrSges(6,1) * t86 - mrSges(6,3) * t278;
t62 = -mrSges(6,2) * t133 - mrSges(6,3) * t244;
t61 = mrSges(5,1) * t133 + mrSges(5,2) * t244;
t60 = pkin(4) * t244 + t312;
t59 = -mrSges(7,2) * t244 + mrSges(7,3) * t133;
t58 = t168 * t330 + t267;
t50 = -t133 * Ifges(5,2) - t227 * Ifges(5,6) + t315;
t49 = -t227 * Ifges(6,4) - Ifges(6,2) * t244 + t130;
t42 = -pkin(4) * t301 + t266;
t40 = -t168 * pkin(5) - t56;
t39 = t133 * pkin(4) + t243;
t38 = t169 * pkin(5) + t241 * qJ(6) + t57;
t35 = t244 * t330 + t312;
t30 = pkin(4) * t113 + t249;
t28 = t86 * mrSges(5,1) - t74;
t27 = -t86 * mrSges(6,2) + t73;
t20 = t133 * t330 + t243;
t13 = -pkin(4) * t299 + t252;
t11 = pkin(4) * t86 + t251;
t10 = qJD(6) * t168 + t113 * t330 + t249;
t9 = -t113 * pkin(5) - t12;
t8 = -t112 * pkin(5) + t241 * qJD(6) + t252 - t269;
t3 = qJD(6) * t133 + t330 * t86 + t251;
t7 = [(Ifges(5,4) * t169 - Ifges(5,2) * t168) * t358 + m(4) * (t116 * t156 + t117 * t157 + t139 * t144 + t140 * t145) + (-t17 * mrSges(7,1) + t32 * mrSges(6,1) + Ifges(6,6) * t344 - Ifges(5,2) * t346 + t384 * t343 + t381 * t345 + t376 / 0.2e1 + t367 * t338 + mrSges(5,1) * t151 - t50 / 0.2e1 + mrSges(7,3) * t20 - mrSges(6,2) * t39 - t37 * mrSges(5,3)) * t113 + (-t14 * mrSges(7,1) - t31 * mrSges(6,1) + Ifges(6,2) * t344 - Ifges(5,4) * t346 - t385 * t343 + t380 * t345 - t375 / 0.2e1 + t368 * t338 - mrSges(5,2) * t151 + mrSges(7,2) * t20 + mrSges(6,3) * t39 + t49 / 0.2e1 - t36 * mrSges(5,3)) * t112 + m(7) * (t1 * t38 + t10 * t20 + t14 * t8 + t17 * t9 + t2 * t40 + t3 * t58) + m(6) * (t11 * t93 + t12 * t32 + t13 * t31 + t30 * t39 + t4 * t56 + t5 * t57) + (t1 * t169 - t168 * t2) * mrSges(7,1) + (t168 * t4 + t169 * t5) * mrSges(6,1) + (-t168 * t6 + t169 * t255) * mrSges(5,3) + m(5) * (t15 * t37 + t151 * t185 + t173 * t196 + t252 * t36 - t255 * t69 + t6 * t70) - t252 * t103 + (t168 * t384 + t169 * t385) * t360 + (t168 * t381 - t169 * t380) * t357 + t196 * t28 + t185 * t61 + t173 * (mrSges(5,1) * t168 + mrSges(5,2) * t169) + t11 * (-mrSges(6,2) * t168 - mrSges(6,3) * t169) + t157 * t171 + t156 * t172 + t3 * (-mrSges(7,2) * t169 + mrSges(7,3) * t168) + t144 * t159 + t145 * t158 - t366 * t241 / 0.2e1 + (((Ifges(4,6) * t334 + t363) * qJD(2) + t364) * qJD(2) + pkin(7) * t163 + t147 * t336 + t148 * t334 + (-t116 * t237 - t117 * t236) * mrSges(4,3) + (-0.2e1 * t356 + (t319 / 0.2e1 - t317 - 0.3e1 / 0.2e1 * Ifges(3,4)) * t240 - t272 * t169 - t271 * t168) * t295) * t240 + t13 * t107 + t15 * t102 + t8 * t104 + t12 * t105 + t9 * t106 + t93 * t27 + t70 * t68 + t38 * t63 + t56 * t64 + t40 * t65 + t57 * t66 + t69 * t67 + t30 * t62 + (-Ifges(6,2) * t169 + Ifges(6,6) * t168) * t359 + t58 * t29 + t10 * t59 + t168 * t394 + t169 * t395 + t168 * t361 + t169 * t362 + (-t116 * mrSges(4,1) + t117 * mrSges(4,2) - Ifges(5,6) * t358 - Ifges(6,4) * t359 + (t182 / 0.2e1 + (t263 * t334 + Ifges(3,5) / 0.2e1) * qJD(2) + t190 * t264 / 0.2e1 + t199 * t265 + t126 * t334 + t125 * t336 + (-t139 * t237 - t140 * t236) * mrSges(4,3) + t270 * pkin(7)) * qJD(2) - t392 * t360 - t393 * t357 + (-0.2e1 * t355 + (0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * t319 + 0.3e1 / 0.2e1 * t317) * t241 + (-0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + Ifges(4,1) * t237 ^ 2 / 0.2e1 + (-0.3e1 / 0.2e1 * t321 + t318) * t236 + (m(4) * pkin(7) + t265) * pkin(7) - t273) * t240) * t295 + t365) * t241; (t147 / 0.2e1 + t117 * mrSges(4,3) + qJD(3) * t158 + qJ(3) * t171) * t237 + (t66 - t67) * t149 + (t11 * t128 + t149 * t5 - t150 * t4 + t372 * t39 + t391 * t32 + (t101 - t42) * t31) * m(6) + (Ifges(5,4) * t161 - Ifges(5,2) * t160 + t178 * t381 - t369 * t380) * t345 + (Ifges(5,4) * t369 - Ifges(5,2) * t178 + t160 * t381 - t161 * t380) * t346 + (t49 - t375) * (-t369 / 0.2e1 + t161 / 0.2e1) + (t178 * t367 - t368 * t369) * t338 + (-Ifges(6,2) * t369 + Ifges(6,6) * t178 + t160 * t384 + t161 * t385) * t344 + (-Ifges(6,2) * t161 + Ifges(6,6) * t160 + t178 * t384 + t369 * t385) * t343 + (t68 - t64) * t150 - m(4) * (t139 * t152 + t140 * t153) + (-t192 * t6 + t193 * t255 - t303 * t36 + t304 * t37) * mrSges(5,3) + t266 * t103 + (t149 * t255 + t150 * t6 - t151 * t184 + t173 * t229 + (-t100 - t55) * t37 + t390 * t36) * m(5) + t229 * t28 + (t192 * t384 + t193 * t385) * t360 + (t192 * t381 - t193 * t380) * t357 + t3 * (-mrSges(7,2) * t193 + mrSges(7,3) * t192) + t173 * (mrSges(5,1) * t192 + mrSges(5,2) * t193) + t11 * (-mrSges(6,2) * t192 - mrSges(6,3) * t193) + t372 * t62 + t373 * t104 + t374 * t106 + (-t50 + t376) * (t178 / 0.2e1 - t160 / 0.2e1) + t377 * t59 + (t1 * t110 + t111 * t2 + t14 * t373 + t17 * t374 + t20 * t377 + t3 * t96) * m(7) - t184 * t61 - pkin(2) * t163 - t152 * t159 + (t160 * t367 - t161 * t368) * t337 - t153 * t158 + t128 * t27 + ((-t182 / 0.2e1 + t313 / 0.2e1 - t230 / 0.2e1 + qJD(1) * t355 + (-t325 / 0.2e1 + t139 * mrSges(4,3) - t199 * mrSges(4,2) - t126 / 0.2e1 + t290 / 0.2e1) * t237 + ((-m(4) * pkin(2) - mrSges(4,1) * t237 - mrSges(3,1)) * qJD(2) - t270) * pkin(7) + (t323 / 0.2e1 + t140 * mrSges(4,3) - t199 * mrSges(4,1) + t125 / 0.2e1 - t289 / 0.2e1 + (t324 / 0.2e1 + pkin(7) * mrSges(4,2) - t322 / 0.2e1) * qJD(2)) * t236) * t241 + ((t263 * t335 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t300 + (pkin(7) * mrSges(3,2) + Ifges(4,5) * t335 - t271 * t192 - t272 * t193 + t363) * qJD(2) + (t356 + (t317 / 0.2e1 + Ifges(3,4) / 0.2e1) * t240) * qJD(1) - t364) * t240) * qJD(1) - t42 * t107 + t110 * t63 + t111 * t65 - t55 * t102 - t41 * t105 + t96 * t29 + t192 * t394 + t193 * t395 + (Ifges(5,4) * t193 - Ifges(5,2) * t192) * t358 + (-Ifges(6,2) * t193 + Ifges(6,6) * t192) * t359 + t192 * t361 + t193 * t362 + (mrSges(6,2) * t304 + mrSges(6,3) * t303) * t39 + (t1 * t193 - t14 * t303 + t17 * t304 - t192 * t2) * mrSges(7,1) + (t148 / 0.2e1 - t116 * mrSges(4,3) - qJD(3) * t159 - qJ(3) * t172) * t236 + (mrSges(7,2) * t303 - mrSges(7,3) * t304) * t20 + (t192 * t4 + t193 * t5 - t303 * t31 - t304 * t32) * mrSges(6,1) + (-mrSges(5,1) * t304 - mrSges(5,2) * t303) * t151 + t305 * t101 + t307 * t100 + m(4) * (t140 * t237 * qJD(3) - t139 * t298 + (-t116 * t236 + t117 * t237) * qJ(3)); -t74 + t73 + t254 * t158 + t190 * t159 + (mrSges(5,1) - mrSges(6,2)) * t86 - (-t102 - t306) * t133 + (-t104 - t305) * t244 + t29 + t163 + (t133 * t17 - t14 * t244 + t3) * m(7) + (-t133 * t32 - t244 * t31 + t11) * m(6) + (t133 * t37 - t244 * t36 + t173) * m(5) + (t139 * t190 + t140 * t254 + t226) * m(4); -t365 + t366 + (-pkin(4) * t5 - qJ(5) * t4 - t31 * t37 + t32 * t371 - t39 * t60) * m(6) + (Ifges(6,2) * t133 + t314 + t50) * t343 + (-t64 + t65) * qJ(5) + (-t133 * t385 + t129 - t315 + t376) * t344 - t330 * t63 + t256 * t106 + t370 * t104 - pkin(4) * t66 - t60 * t62 - t35 * t59 + t306 * qJD(5) + (-t305 + t326) * t37 + (-t307 + t327) * t36 + (t133 * t14 + t17 * t244) * mrSges(7,1) + (t133 * t31 - t244 * t32) * mrSges(6,1) + (-Ifges(5,2) * t244 - t131 + t375) * t345 + (t244 * t381 + t130 - t316 + t49) * t346 - t151 * (mrSges(5,1) * t244 - mrSges(5,2) * t133) - t20 * (mrSges(7,2) * t133 + mrSges(7,3) * t244) - t39 * (-mrSges(6,2) * t244 + mrSges(6,3) * t133) + (t133 * t368 + t244 * t367) * t337 + (qJ(5) * t2 - t1 * t330 + t370 * t14 + t17 * t388 - t20 * t35) * m(7); t306 * t227 + (t59 + t62) * t244 + t63 + t66 + (t17 * t227 + t20 * t244 + t1) * m(7) + (-t227 * t32 + t244 * t39 + t5) * m(6); -t227 * t104 - t133 * t59 + 0.2e1 * (t2 / 0.2e1 + t20 * t346 + t14 * t338) * m(7) + t65;];
tauc  = t7(:);
