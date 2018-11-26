% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2018-11-23 15:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:29:42
% EndTime: 2018-11-23 15:29:58
% DurationCPUTime: 16.77s
% Computational Cost: add. (7954->675), mult. (19898->931), div. (0->0), fcn. (14068->10), ass. (0->276)
t372 = Ifges(6,1) + Ifges(7,1);
t371 = Ifges(7,4) + Ifges(6,5);
t370 = Ifges(7,5) - Ifges(6,4);
t223 = sin(qJ(5));
t224 = sin(qJ(4));
t227 = cos(qJ(5));
t228 = cos(qJ(4));
t183 = t223 * t224 - t227 * t228;
t357 = qJD(4) + qJD(5);
t126 = t357 * t183;
t229 = cos(qJ(3));
t237 = t183 * t229;
t161 = qJD(2) * t237;
t383 = t126 - t161;
t184 = t223 * t228 + t224 * t227;
t127 = t357 * t184;
t293 = qJD(2) * t229;
t160 = t184 * t293;
t300 = t127 - t160;
t369 = Ifges(7,2) + Ifges(6,3);
t368 = -Ifges(6,6) + Ifges(7,6);
t347 = -pkin(10) - pkin(9);
t198 = t347 * t224;
t199 = t347 * t228;
t146 = t198 * t223 - t199 * t227;
t278 = qJD(4) * t347;
t188 = t224 * t278;
t260 = t228 * t278;
t226 = sin(qJ(2));
t221 = sin(pkin(6));
t297 = qJD(1) * t221;
t277 = t226 * t297;
t191 = qJD(2) * pkin(8) + t277;
t225 = sin(qJ(3));
t222 = cos(pkin(6));
t296 = qJD(1) * t222;
t154 = -t225 * t191 + t229 * t296;
t256 = pkin(3) * t225 - pkin(9) * t229;
t187 = t256 * qJD(2);
t102 = -t154 * t224 + t228 * t187;
t304 = t228 * t229;
t243 = pkin(4) * t225 - pkin(10) * t304;
t87 = qJD(2) * t243 + t102;
t103 = t228 * t154 + t224 * t187;
t273 = t224 * t293;
t95 = -pkin(10) * t273 + t103;
t379 = -qJD(5) * t146 + (t260 - t87) * t227 + (-t188 + t95) * t223;
t292 = qJD(3) * t224;
t294 = qJD(2) * t225;
t182 = t228 * t294 + t292;
t290 = qJD(3) * t228;
t239 = t224 * t294 - t290;
t119 = t223 * t182 + t227 * t239;
t114 = Ifges(6,4) * t119;
t210 = qJD(4) - t293;
t200 = qJD(5) + t210;
t234 = t227 * t182 - t223 * t239;
t318 = Ifges(7,5) * t119;
t367 = t371 * t200 + t234 * t372 - t114 + t318;
t382 = t367 / 0.2e1;
t283 = qJD(2) * qJD(3);
t266 = t225 * t283;
t282 = qJD(3) * qJD(4);
t287 = qJD(4) * t225;
t289 = qJD(3) * t229;
t360 = -t224 * t287 + t228 * t289;
t147 = qJD(2) * t360 + t228 * t282;
t272 = t224 * t289;
t286 = qJD(4) * t228;
t376 = t225 * t286 + t272;
t148 = -qJD(2) * t376 - t224 * t282;
t47 = -qJD(5) * t119 + t227 * t147 + t223 * t148;
t48 = qJD(5) * t234 + t223 * t147 - t227 * t148;
t381 = t371 * t266 + t370 * t48 + t372 * t47;
t155 = t229 * t191 + t225 * t296;
t124 = pkin(4) * t273 + t155;
t288 = qJD(4) * t224;
t380 = pkin(4) * t288 + pkin(5) * t300 + qJ(6) * t383 - qJD(6) * t184 - t124;
t378 = pkin(5) * t294 - t379;
t230 = cos(qJ(2));
t303 = t229 * t230;
t149 = (-t224 * t303 + t226 * t228) * t297;
t150 = (t224 * t226 + t228 * t303) * t297;
t193 = -pkin(3) * t229 - pkin(9) * t225 - pkin(2);
t181 = t228 * t193;
t305 = t225 * t228;
t329 = pkin(8) * t224;
t123 = -pkin(10) * t305 + t181 + (-pkin(4) - t329) * t229;
t212 = pkin(8) * t304;
t159 = t224 * t193 + t212;
t306 = t224 * t225;
t133 = -pkin(10) * t306 + t159;
t361 = t223 * t123 + t227 * t133;
t189 = t256 * qJD(3);
t291 = qJD(3) * t225;
t298 = t228 * t189 + t291 * t329;
t72 = t243 * qJD(3) + (-t212 + (pkin(10) * t225 - t193) * t224) * qJD(4) + t298;
t96 = t224 * t189 + t193 * t286 + (-t225 * t290 - t229 * t288) * pkin(8);
t78 = -pkin(10) * t376 + t96;
t365 = -qJD(5) * t361 + (-t149 + t72) * t227 + (t150 - t78) * t223;
t350 = t48 / 0.2e1;
t377 = Ifges(7,3) * t350;
t69 = pkin(5) * t234 + qJ(6) * t119;
t231 = qJD(2) ^ 2;
t352 = t47 / 0.2e1;
t341 = t147 / 0.2e1;
t340 = t148 / 0.2e1;
t374 = t266 / 0.2e1;
t366 = -pkin(5) * t291 - t365;
t364 = t96 - t150;
t363 = -qJD(4) * t159 - t149 + t298;
t162 = t184 * t225;
t141 = qJD(3) * pkin(9) + t155;
t276 = t230 * t297;
t157 = qJD(2) * t193 - t276;
t89 = t228 * t141 + t224 * t157;
t68 = -pkin(10) * t239 + t89;
t315 = t223 * t68;
t88 = -t141 * t224 + t228 * t157;
t67 = -pkin(10) * t182 + t88;
t58 = pkin(4) * t210 + t67;
t25 = t227 * t58 - t315;
t362 = qJD(6) - t25;
t359 = t369 * t266 + t368 * t48 + t371 * t47;
t295 = qJD(2) * t221;
t274 = t230 * t295;
t235 = qJD(1) * (qJD(3) * t222 + t274);
t109 = -t191 * t291 + t229 * t235;
t151 = (t189 + t277) * qJD(2);
t31 = t228 * t109 - t141 * t288 + t224 * t151 + t157 * t286;
t32 = -t89 * qJD(4) - t109 * t224 + t228 * t151;
t358 = -t224 * t32 + t228 * t31;
t322 = Ifges(5,4) * t182;
t107 = -Ifges(5,2) * t239 + Ifges(5,6) * t210 + t322;
t179 = Ifges(5,4) * t239;
t108 = t182 * Ifges(5,1) + t210 * Ifges(5,5) - t179;
t140 = -qJD(3) * pkin(3) - t154;
t248 = t224 * t89 + t228 * t88;
t317 = Ifges(5,6) * t224;
t251 = Ifges(5,5) * t228 - t317;
t321 = Ifges(5,4) * t224;
t255 = Ifges(5,1) * t228 - t321;
t332 = t228 / 0.2e1;
t333 = t210 / 0.2e1;
t337 = t182 / 0.2e1;
t356 = -t248 * mrSges(5,3) + t140 * (mrSges(5,1) * t224 + mrSges(5,2) * t228) + t255 * t337 + t251 * t333 - t224 * t107 / 0.2e1 + t108 * t332;
t24 = pkin(4) * t266 - pkin(10) * t147 + t32;
t313 = t227 * t68;
t26 = t223 * t58 + t313;
t28 = pkin(10) * t148 + t31;
t6 = -qJD(5) * t26 - t223 * t28 + t227 * t24;
t105 = pkin(4) * t239 + t140;
t22 = qJ(6) * t200 + t26;
t37 = t119 * pkin(5) - qJ(6) * t234 + t105;
t113 = Ifges(7,5) * t234;
t52 = Ifges(7,6) * t200 + Ifges(7,3) * t119 + t113;
t319 = Ifges(6,4) * t234;
t55 = -Ifges(6,2) * t119 + Ifges(6,6) * t200 + t319;
t355 = t105 * mrSges(6,1) + t37 * mrSges(7,1) + t52 / 0.2e1 - t55 / 0.2e1 - t22 * mrSges(7,2) - t26 * mrSges(6,3);
t354 = Ifges(7,5) * t352 + Ifges(7,6) * t374 + t377;
t353 = -t47 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t350 - Ifges(6,6) * t266 / 0.2e1;
t351 = -t48 / 0.2e1;
t348 = Ifges(5,1) * t341 + Ifges(5,4) * t340 + Ifges(5,5) * t374;
t346 = -t119 / 0.2e1;
t345 = t119 / 0.2e1;
t343 = -t234 / 0.2e1;
t342 = t234 / 0.2e1;
t338 = -t182 / 0.2e1;
t335 = -t200 / 0.2e1;
t334 = t200 / 0.2e1;
t330 = m(6) * t105;
t328 = qJD(3) / 0.2e1;
t38 = -mrSges(7,2) * t48 + mrSges(7,3) * t266;
t41 = -mrSges(6,2) * t266 - mrSges(6,3) * t48;
t327 = t38 + t41;
t39 = mrSges(6,1) * t266 - mrSges(6,3) * t47;
t40 = -mrSges(7,1) * t266 + t47 * mrSges(7,2);
t326 = t40 - t39;
t36 = t223 * t87 + t227 * t95;
t98 = -mrSges(7,2) * t119 + mrSges(7,3) * t200;
t324 = mrSges(6,3) * t119;
t99 = -mrSges(6,2) * t200 - t324;
t325 = t98 + t99;
t323 = mrSges(6,3) * t234;
t320 = Ifges(5,4) * t228;
t316 = Ifges(5,3) * t210;
t311 = Ifges(4,6) * qJD(3);
t110 = t191 * t289 + t225 * t235;
t308 = t221 * t226;
t166 = -t222 * t229 + t225 * t308;
t310 = t110 * t166;
t309 = t110 * t225;
t307 = t221 * t230;
t100 = mrSges(6,1) * t200 - t323;
t101 = -mrSges(7,1) * t200 + mrSges(7,2) * t234;
t302 = t100 - t101;
t299 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t239 + t182 * mrSges(5,2) + mrSges(4,3) * t294;
t190 = pkin(4) * t306 + t225 * pkin(8);
t285 = qJD(5) * t223;
t284 = qJD(5) * t227;
t279 = Ifges(5,5) * t147 + Ifges(5,6) * t148 + Ifges(5,3) * t266;
t156 = pkin(4) * t376 + pkin(8) * t289;
t217 = -pkin(4) * t228 - pkin(3);
t275 = t226 * t295;
t265 = t229 * t283;
t263 = t289 / 0.2e1;
t262 = -t287 / 0.2e1;
t70 = mrSges(7,1) * t119 - mrSges(7,3) * t234;
t71 = mrSges(6,1) * t119 + mrSges(6,2) * t234;
t261 = t70 + t71 + t299;
t259 = t225 * t276;
t254 = Ifges(5,1) * t224 + t320;
t253 = -Ifges(5,2) * t224 + t320;
t252 = Ifges(5,2) * t228 + t321;
t247 = t109 * t229 + t309;
t73 = t123 * t227 - t133 * t223;
t167 = t222 * t225 + t229 * t308;
t131 = -t167 * t224 - t228 * t307;
t242 = -t167 * t228 + t224 * t307;
t245 = t227 * t131 + t223 * t242;
t77 = t131 * t223 - t227 * t242;
t244 = t227 * t198 + t199 * t223;
t5 = t223 * t24 + t227 * t28 + t58 * t284 - t285 * t68;
t241 = t224 * t253;
t240 = t228 * t253;
t12 = t123 * t284 - t133 * t285 + t223 * t72 + t227 * t78;
t238 = (mrSges(4,1) * t225 + mrSges(4,2) * t229) * qJD(3);
t2 = qJ(6) * t266 + qJD(6) * t200 + t5;
t3 = -pkin(5) * t266 - t6;
t233 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) + t359;
t79 = -pkin(4) * t148 + t110;
t218 = Ifges(4,4) * t293;
t216 = -pkin(4) * t227 - pkin(5);
t214 = pkin(4) * t223 + qJ(6);
t205 = pkin(4) * t284 + qJD(6);
t197 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t293;
t192 = -qJD(2) * pkin(2) - t276;
t175 = qJD(2) * t238;
t169 = Ifges(4,1) * t294 + Ifges(4,5) * qJD(3) + t218;
t168 = t311 + (Ifges(4,4) * t225 + t229 * Ifges(4,2)) * qJD(2);
t163 = t183 * t225;
t158 = -t229 * t329 + t181;
t153 = mrSges(5,1) * t210 - mrSges(5,3) * t182;
t152 = -t210 * mrSges(5,2) - mrSges(5,3) * t239;
t130 = qJD(3) * t167 + t225 * t274;
t129 = -qJD(3) * t166 + t229 * t274;
t117 = -mrSges(5,2) * t266 + mrSges(5,3) * t148;
t116 = mrSges(5,1) * t266 - mrSges(5,3) * t147;
t115 = pkin(5) * t183 - qJ(6) * t184 + t217;
t106 = Ifges(5,5) * t182 - Ifges(5,6) * t239 + t316;
t94 = pkin(5) * t162 + qJ(6) * t163 + t190;
t93 = -mrSges(5,1) * t148 + mrSges(5,2) * t147;
t91 = t149 * t223 + t150 * t227;
t85 = qJD(5) * t244 + t227 * t188 + t223 * t260;
t82 = t147 * Ifges(5,4) + t148 * Ifges(5,2) + Ifges(5,6) * t266;
t81 = -t285 * t306 + (t305 * t357 + t272) * t227 + t360 * t223;
t80 = -qJD(3) * t237 - t162 * t357;
t66 = pkin(5) * t229 - t73;
t65 = -qJ(6) * t229 + t361;
t62 = qJD(4) * t242 - t129 * t224 + t228 * t275;
t61 = qJD(4) * t131 + t129 * t228 + t224 * t275;
t54 = Ifges(7,4) * t234 + t200 * Ifges(7,2) + t119 * Ifges(7,6);
t53 = Ifges(6,5) * t234 - t119 * Ifges(6,6) + t200 * Ifges(6,3);
t50 = pkin(4) * t182 + t69;
t33 = qJ(6) * t294 + t36;
t30 = t227 * t67 - t315;
t29 = t223 * t67 + t313;
t21 = pkin(5) * t81 - qJ(6) * t80 + qJD(6) * t163 + t156;
t20 = -pkin(5) * t200 + t362;
t19 = mrSges(6,1) * t48 + mrSges(6,2) * t47;
t18 = mrSges(7,1) * t48 - mrSges(7,3) * t47;
t10 = qJ(6) * t291 - qJD(6) * t229 + t12;
t9 = qJD(5) * t77 + t223 * t61 - t227 * t62;
t8 = qJD(5) * t245 + t223 * t62 + t227 * t61;
t7 = pkin(5) * t48 - qJ(6) * t47 - qJD(6) * t234 + t79;
t1 = [-t167 * mrSges(4,3) * t266 + t131 * t116 - t242 * t117 + t129 * t197 + t61 * t152 + t62 * t153 - t302 * t9 + t325 * t8 + t327 * t77 - t326 * t245 + ((-mrSges(3,2) * t231 - t175) * t230 + (-mrSges(4,1) * t229 + mrSges(4,2) * t225 - mrSges(3,1)) * t231 * t226) * t221 + (mrSges(4,3) * t265 + t18 + t19 + t93) * t166 + t261 * t130 + m(5) * (t130 * t140 + t131 * t32 - t242 * t31 + t61 * t89 + t62 * t88 + t310) + m(6) * (t105 * t130 + t166 * t79 + t245 * t6 - t25 * t9 + t26 * t8 + t5 * t77) + m(7) * (t130 * t37 + t166 * t7 + t2 * t77 + t20 * t9 + t22 * t8 - t245 * t3) + m(4) * (t109 * t167 + t310 + t129 * t155 - t130 * t154 + (t192 - t276) * t275); (-(t192 * t226 + (-t154 * t225 + t155 * t229) * t230) * t297 - pkin(2) * qJD(1) * t275) * m(4) + (t106 + t54 + t53 + (t225 * t251 - t371 * t163 + t368 * t162 + (-Ifges(5,3) - t369) * t229) * qJD(2)) * t291 / 0.2e1 + t361 * t41 + (t190 * t79 + t5 * t361 + t6 * t73 + (t12 - t91) * t26 + t365 * t25 + (t156 - t259) * t105) * m(6) + (Ifges(6,4) * t80 + Ifges(6,6) * t291) * t346 + (-Ifges(7,5) * t163 - Ifges(7,6) * t229) * t350 + (-Ifges(6,4) * t163 - Ifges(6,6) * t229) * t351 + (t105 * t80 - t163 * t79 + t229 * t5 - t26 * t291) * mrSges(6,2) + (t163 * t7 - t2 * t229 + t22 * t291 - t37 * t80) * mrSges(7,3) + t6 * (-mrSges(6,1) * t229 + mrSges(6,3) * t163) + t3 * (mrSges(7,1) * t229 - mrSges(7,2) * t163) + (t224 * t262 + t228 * t263) * t108 + (-t272 / 0.2e1 + t228 * t262) * t107 + t94 * t18 + t10 * t98 + t12 * t99 + t73 * t39 + t65 * t38 + t66 * t40 + t21 * t70 - (t279 + t359) * t229 / 0.2e1 + ((-t154 * t229 - t155 * t225) * qJD(3) + t247) * mrSges(4,3) - t82 * t306 / 0.2e1 - t325 * t91 + (-t229 * t197 - t225 * t261) * t276 + t363 * t153 + (-t140 * t259 + t158 * t32 + t159 * t31 + t363 * t88 + t364 * t89) * m(5) + t364 * t152 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t225 * t265 + ((t140 * t290 + t31) * t229 + (-qJD(3) * t89 + t110 * t228 - t140 * t288) * t225) * mrSges(5,2) - t168 * t291 / 0.2e1 + t20 * (-mrSges(7,1) * t291 + mrSges(7,2) * t80) + t25 * (mrSges(6,1) * t291 - mrSges(6,3) * t80) + ((t140 * t292 - t32) * t229 + (qJD(3) * t88 + t110 * t224 + t140 * t286) * t225) * mrSges(5,1) - t239 * (-t252 * t287 + (Ifges(5,6) * t225 + t229 * t253) * qJD(3)) / 0.2e1 + (-t248 * t289 + (-t31 * t224 - t32 * t228 + (t224 * t88 - t228 * t89) * qJD(4)) * t225) * mrSges(5,3) + t169 * t263 + qJD(3) ^ 2 * (Ifges(4,5) * t229 - Ifges(4,6) * t225) / 0.2e1 + (m(5) * (t140 * t289 + t309) + m(4) * (-t154 * t289 - t155 * t291 + t247) + t225 * t93 + (-t197 * t225 + t229 * t299) * qJD(3)) * pkin(8) + t365 * t100 + (t2 * t65 + t3 * t66 + t7 * t94 + (t21 - t259) * t37 + (t10 - t91) * t22 + t366 * t20) * m(7) + t366 * t101 + (-Ifges(6,2) * t346 + Ifges(7,3) * t345 + t334 * t368 + t342 * t370 + t355) * t81 + (t291 * t369 + t371 * t80) * t334 + (t291 * t371 + t372 * t80) * t342 + (-t163 * t372 - t229 * t371) * t352 + t156 * t71 + (t79 * mrSges(6,1) + t7 * mrSges(7,1) - t2 * mrSges(7,2) - t5 * mrSges(6,3) - Ifges(6,2) * t351 + t352 * t370 + t353 + t354 + t377) * t162 + t158 * t116 + t159 * t117 - t381 * t163 / 0.2e1 - pkin(2) * t175 + (0.3e1 / 0.2e1 * t229 ^ 2 - 0.3e1 / 0.2e1 * t225 ^ 2) * Ifges(4,4) * t283 + t190 * t19 + (Ifges(7,5) * t80 + Ifges(7,6) * t291) * t345 + t80 * t382 + t192 * t238 + ((-Ifges(5,5) * t224 - Ifges(5,6) * t228) * t287 + (Ifges(5,3) * t225 + t229 * t251) * qJD(3)) * t333 + (-t254 * t287 + (Ifges(5,5) * t225 + t229 * t255) * qJD(3)) * t337 + (-Ifges(5,6) * t229 + t225 * t253) * t340 + (-Ifges(5,5) * t229 + t225 * t255) * t341 + t305 * t348; (-Ifges(6,4) * t161 - Ifges(7,5) * t126 - Ifges(6,2) * t160 + Ifges(7,3) * t127) * t345 + (-Ifges(6,4) * t126 - Ifges(7,5) * t161 - Ifges(6,2) * t127 + Ifges(7,3) * t160) * t346 + (t52 - t55) * (t127 / 0.2e1 - t160 / 0.2e1) - t326 * t244 + (-mrSges(5,1) * t228 + mrSges(5,2) * t224 - mrSges(4,1)) * t110 + t115 * t18 - t109 * mrSges(4,2) - pkin(3) * t93 - t33 * t98 - t36 * t99 + t358 * mrSges(5,3) + m(5) * (-pkin(3) * t110 + pkin(9) * t358) + t325 * t85 + t327 * t146 + (-t116 * t224 + t117 * t228) * pkin(9) - t299 * t155 + (t240 * t328 + (t71 + t330) * t224 * pkin(4) + (-m(5) * t248 - t224 * t152 - t228 * t153) * pkin(9) + t356) * qJD(4) - m(5) * (t102 * t88 + t103 * t89 + t140 * t155) - t124 * t71 + ((t154 * mrSges(4,3) - t218 / 0.2e1 - t192 * mrSges(4,2) - t169 / 0.2e1 + (Ifges(4,5) / 0.2e1 - t240 / 0.2e1) * qJD(3) - t356) * t229 + ((Ifges(4,4) / 0.2e1 + t317 / 0.2e1) * t294 + t155 * mrSges(4,3) - t25 * mrSges(6,1) + t20 * mrSges(7,1) - t22 * mrSges(7,3) + t26 * mrSges(6,2) - t192 * mrSges(4,1) - t88 * mrSges(5,1) + t89 * mrSges(5,2) - t54 / 0.2e1 - t106 / 0.2e1 + t168 / 0.2e1 - t53 / 0.2e1 - t311 / 0.2e1 - t316 / 0.2e1 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t200 + (-Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1) * t234 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t119 - qJD(4) * t241 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1 + t241 / 0.2e1) * t293 + (t292 / 0.2e1 + t338) * Ifges(5,5) + (t183 * t368 + t184 * t371) * t328) * t225) * qJD(2) - t103 * t152 - t102 * t153 + t378 * t101 + (-t105 * t124 + t244 * t6 + t146 * t5 + t217 * t79 + (-t36 + t85) * t26 + t379 * t25) * m(6) + t379 * t100 + t380 * t70 + (t115 * t7 - t244 * t3 + t146 * t2 + t380 * t37 + (-t33 + t85) * t22 + t378 * t20) * m(7) + t367 * (-t126 / 0.2e1 + t161 / 0.2e1) + t381 * t184 / 0.2e1 + (-t126 * t371 + t127 * t368) * t334 + (t160 * t368 - t161 * t371) * t335 + (-t126 * t372 + t127 * t370) * t342 + (t160 * t370 - t161 * t372) * t343 + (t183 * t370 + t184 * t372) * t352 + t7 * (mrSges(7,1) * t183 - mrSges(7,3) * t184) + t79 * (mrSges(6,1) * t183 + mrSges(6,2) * t184) - t154 * t197 + t217 * t19 + t82 * t332 + t252 * t340 + t254 * t341 + t224 * t348 + (Ifges(7,5) * t184 + Ifges(7,3) * t183) * t350 + (Ifges(6,4) * t184 - Ifges(6,2) * t183) * t351 + t183 * t353 + t183 * t354 + (mrSges(6,1) * t300 - mrSges(6,2) * t383) * t105 + (-t183 * t2 + t184 * t3 - t20 * t383 - t22 * t300) * mrSges(7,2) + (mrSges(7,1) * t300 + mrSges(7,3) * t383) * t37 + (-t183 * t5 - t184 * t6 + t25 * t383 - t26 * t300) * mrSges(6,3); (t2 * t214 - t20 * t29 + t216 * t3 - t37 * t50 + (t205 - t30) * t22) * m(7) - t50 * t70 - t31 * mrSges(5,2) + t32 * mrSges(5,1) - t325 * t30 + t302 * t29 - m(6) * (-t25 * t29 + t26 * t30) + t279 + t233 + (t89 * t182 - t239 * t88) * mrSges(5,3) - t140 * (t182 * mrSges(5,1) - mrSges(5,2) * t239) - t210 * (-Ifges(5,5) * t239 - Ifges(5,6) * t182) / 0.2e1 + (-Ifges(6,2) * t345 + Ifges(7,3) * t346 + t335 * t368 + t343 * t370 - t355) * t234 - t88 * t152 + t89 * t153 + t205 * t98 + t214 * t38 + t216 * t40 + (m(6) * (t223 * t5 + t227 * t6 - t25 * t285 + t26 * t284) + m(7) * t20 * t285 + 0.2e1 * t330 * t338 - t182 * t71 + t223 * t41 + t227 * t39 + (-t223 * t302 + t227 * t99) * qJD(5)) * pkin(4) + t107 * t337 + (-Ifges(5,1) * t239 - t322) * t338 + (-Ifges(5,2) * t182 + t108 - t179) * t239 / 0.2e1 + (t105 * mrSges(6,2) + t20 * mrSges(7,2) - t25 * mrSges(6,3) - t37 * mrSges(7,3) - Ifges(6,4) * t345 - Ifges(7,5) * t346 - t371 * t335 - t372 * t343 + t382) * t119; qJD(6) * t98 - t69 * t70 + qJ(6) * t38 - pkin(5) * t40 + (t302 + t323) * t26 + (-t324 - t325) * t25 + t233 + (t119 * t20 + t22 * t234) * mrSges(7,2) - t37 * (mrSges(7,1) * t234 + mrSges(7,3) * t119) + (Ifges(7,3) * t234 - t318) * t346 + t55 * t342 - t105 * (mrSges(6,1) * t234 - mrSges(6,2) * t119) + (-t119 * t371 + t234 * t368) * t335 + (-pkin(5) * t3 + qJ(6) * t2 - t20 * t26 + t22 * t362 - t37 * t69) * m(7) + (-Ifges(6,2) * t234 - t114 + t367) * t345 + (-t119 * t372 + t113 - t319 + t52) * t343; t234 * t70 - t200 * t98 + 0.2e1 * (t3 / 0.2e1 + t37 * t342 + t22 * t335) * m(7) + t40;];
tauc  = t1(:);
