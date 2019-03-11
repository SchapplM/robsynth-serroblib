% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:15
% EndTime: 2019-03-09 05:00:36
% DurationCPUTime: 12.30s
% Computational Cost: add. (9942->631), mult. (23894->882), div. (0->0), fcn. (16263->10), ass. (0->285)
t246 = sin(pkin(11));
t248 = cos(pkin(11));
t251 = sin(qJ(4));
t254 = cos(qJ(4));
t213 = t246 * t254 + t248 * t251;
t255 = cos(qJ(3));
t261 = t213 * t255;
t174 = qJD(1) * t261;
t200 = t213 * qJD(4);
t399 = t174 - t200;
t263 = t246 * t251 - t248 * t254;
t260 = t263 * t255;
t175 = qJD(1) * t260;
t201 = t263 * qJD(4);
t303 = -t175 + t201;
t330 = -qJ(5) - pkin(8);
t278 = qJD(4) * t330;
t291 = qJD(5) * t254;
t196 = t251 * t278 + t291;
t197 = -qJD(5) * t251 + t254 * t278;
t133 = -t196 * t246 + t248 * t197;
t241 = sin(pkin(10)) * pkin(1) + pkin(7);
t225 = t241 * qJD(1);
t252 = sin(qJ(3));
t188 = qJD(2) * t255 - t252 * t225;
t276 = pkin(3) * t252 - pkin(8) * t255;
t222 = t276 * qJD(1);
t142 = -t251 * t188 + t254 * t222;
t306 = t254 * t255;
t262 = pkin(4) * t252 - qJ(5) * t306;
t114 = qJD(1) * t262 + t142;
t143 = t254 * t188 + t251 * t222;
t299 = qJD(1) * t255;
t285 = t251 * t299;
t123 = -qJ(5) * t285 + t143;
t60 = t248 * t114 - t123 * t246;
t398 = -t60 + t133;
t134 = t248 * t196 + t246 * t197;
t61 = t246 * t114 + t248 * t123;
t397 = -t61 + t134;
t250 = sin(qJ(6));
t253 = cos(qJ(6));
t239 = qJD(4) - t299;
t296 = qJD(3) * t254;
t300 = qJD(1) * t252;
t219 = -t251 * t300 + t296;
t220 = qJD(3) * t251 + t254 * t300;
t153 = t219 * t246 + t220 * t248;
t385 = pkin(9) * t153;
t189 = t252 * qJD(2) + t255 * t225;
t173 = qJD(3) * pkin(8) + t189;
t286 = -cos(pkin(10)) * pkin(1) - pkin(2);
t209 = -pkin(3) * t255 - t252 * pkin(8) + t286;
t178 = t209 * qJD(1);
t120 = -t173 * t251 + t254 * t178;
t98 = -qJ(5) * t220 + t120;
t88 = pkin(4) * t239 + t98;
t121 = t173 * t254 + t178 * t251;
t99 = qJ(5) * t219 + t121;
t92 = t246 * t99;
t44 = t248 * t88 - t92;
t32 = pkin(5) * t239 - t385 + t44;
t277 = t248 * t219 - t220 * t246;
t374 = pkin(9) * t277;
t316 = t248 * t99;
t45 = t246 * t88 + t316;
t35 = t45 + t374;
t10 = t250 * t32 + t253 * t35;
t297 = qJD(3) * t252;
t279 = qJD(1) * t297;
t236 = Ifges(7,3) * t279;
t87 = t153 * t253 + t250 * t277;
t336 = Ifges(7,4) * t87;
t232 = qJD(6) + t239;
t342 = -t232 / 0.2e1;
t360 = -t87 / 0.2e1;
t377 = -t153 * t250 + t253 * t277;
t362 = -t377 / 0.2e1;
t81 = Ifges(7,4) * t377;
t40 = Ifges(7,1) * t87 + Ifges(7,5) * t232 + t81;
t172 = -qJD(3) * pkin(3) - t188;
t144 = -t219 * pkin(4) + qJD(5) + t172;
t89 = -pkin(5) * t277 + t144;
t9 = -t250 * t35 + t253 * t32;
t396 = t236 + (Ifges(7,5) * t377 - Ifges(7,6) * t87) * t342 + (t10 * t87 + t377 * t9) * mrSges(7,3) + (-Ifges(7,2) * t87 + t40 + t81) * t362 - t89 * (mrSges(7,1) * t87 + mrSges(7,2) * t377) + (Ifges(7,1) * t377 - t336) * t360;
t395 = -pkin(5) * t300 + pkin(9) * t303 + t398;
t394 = pkin(9) * t399 + t397;
t392 = -Ifges(4,1) / 0.2e1;
t39 = Ifges(7,2) * t377 + Ifges(7,6) * t232 + t336;
t391 = t39 / 0.2e1;
t244 = Ifges(4,4) * t299;
t390 = -t244 / 0.2e1;
t282 = Ifges(4,5) * qJD(3) / 0.2e1;
t229 = t330 * t251;
t230 = t330 * t254;
t158 = t248 * t229 + t230 * t246;
t135 = -pkin(9) * t213 + t158;
t159 = t246 * t229 - t248 * t230;
t136 = -pkin(9) * t263 + t159;
t69 = t135 * t250 + t136 * t253;
t384 = -qJD(6) * t69 - t250 * t394 + t253 * t395;
t68 = t135 * t253 - t136 * t250;
t383 = qJD(6) * t68 + t250 * t395 + t253 * t394;
t382 = Ifges(6,4) * t153;
t160 = pkin(4) * t285 + t189;
t294 = qJD(4) * t251;
t381 = pkin(4) * t294 - pkin(5) * t399 - t160;
t320 = t220 * Ifges(5,4);
t139 = t219 * Ifges(5,2) + t239 * Ifges(5,6) + t320;
t210 = Ifges(5,4) * t219;
t140 = Ifges(5,1) * t220 + Ifges(5,5) * t239 + t210;
t267 = t120 * t254 + t121 * t251;
t327 = Ifges(5,4) * t254;
t271 = -Ifges(5,2) * t251 + t327;
t328 = Ifges(5,4) * t251;
t273 = Ifges(5,1) * t254 - t328;
t274 = mrSges(5,1) * t251 + mrSges(5,2) * t254;
t325 = Ifges(5,6) * t251;
t326 = Ifges(5,5) * t254;
t337 = t254 / 0.2e1;
t338 = -t251 / 0.2e1;
t339 = t239 / 0.2e1;
t343 = t220 / 0.2e1;
t256 = -t267 * mrSges(5,3) + t172 * t274 + t219 * t271 / 0.2e1 + t273 * t343 + (-t325 + t326) * t339 + t139 * t338 + t140 * t337;
t380 = t188 * mrSges(4,3) + t300 * t392 - t256 - t282 + t390;
t292 = qJD(4) * t254;
t295 = qJD(3) * t255;
t379 = t251 * t295 + t252 * t292;
t290 = qJD(3) * qJD(4);
t293 = qJD(4) * t252;
t165 = t254 * t290 + (-t251 * t293 + t254 * t295) * qJD(1);
t166 = -qJD(1) * t379 - t251 * t290;
t110 = t165 * t248 + t166 * t246;
t176 = t188 * qJD(3);
t223 = t276 * qJD(3);
t208 = qJD(1) * t223;
t66 = -qJD(4) * t121 - t176 * t251 + t254 * t208;
t46 = pkin(4) * t279 - qJ(5) * t165 - qJD(5) * t220 + t66;
t65 = -t173 * t294 + t254 * t176 + t178 * t292 + t251 * t208;
t52 = qJ(5) * t166 + qJD(5) * t219 + t65;
t15 = -t246 * t52 + t248 * t46;
t11 = pkin(5) * t279 - pkin(9) * t110 + t15;
t109 = -t165 * t246 + t166 * t248;
t16 = t246 * t46 + t248 * t52;
t14 = pkin(9) * t109 + t16;
t2 = qJD(6) * t9 + t11 * t250 + t14 * t253;
t3 = -qJD(6) * t10 + t11 * t253 - t14 * t250;
t30 = qJD(6) * t377 + t109 * t250 + t110 * t253;
t31 = -qJD(6) * t87 + t109 * t253 - t110 * t250;
t378 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t30 + Ifges(7,6) * t31;
t77 = Ifges(6,2) * t277 + Ifges(6,6) * t239 + t382;
t376 = t77 / 0.2e1;
t353 = -t277 / 0.2e1;
t281 = -Ifges(4,6) * qJD(3) / 0.2e1;
t57 = -t109 * mrSges(6,1) + t110 * mrSges(6,2);
t8 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t375 = -t57 - t8;
t373 = Ifges(6,4) * t277;
t242 = pkin(4) * t248 + pkin(5);
t334 = pkin(4) * t246;
t195 = t242 * t250 + t253 * t334;
t50 = -t246 * t98 - t316;
t36 = t50 - t374;
t51 = t248 * t98 - t92;
t37 = t51 - t385;
t372 = -t195 * qJD(6) + t250 * t37 - t253 * t36;
t194 = t242 * t253 - t250 * t334;
t371 = t194 * qJD(6) - t250 * t36 - t253 * t37;
t221 = t241 * t306;
t155 = t251 * t209 + t221;
t268 = -t251 * t66 + t254 * t65;
t227 = t286 * qJD(1);
t289 = Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t329 = Ifges(4,4) * t252;
t366 = -t289 * t239 - t120 * mrSges(5,1) - t227 * mrSges(4,1) - t44 * mrSges(6,1) - t9 * mrSges(7,1) - t220 * Ifges(5,5) - t219 * Ifges(5,6) - t281 + (Ifges(4,2) * t255 + t329) * qJD(1) / 0.2e1 - t232 * Ifges(7,3) - t87 * Ifges(7,5) - t377 * Ifges(7,6) - t153 * Ifges(6,5) - t277 * Ifges(6,6) + t10 * mrSges(7,2) + t121 * mrSges(5,2) + t45 * mrSges(6,2) - (Ifges(5,3) + Ifges(6,3)) * t339;
t365 = -t66 * mrSges(5,1) - t15 * mrSges(6,1) + t65 * mrSges(5,2) + t16 * mrSges(6,2) - Ifges(5,5) * t165 - Ifges(6,5) * t110 - Ifges(5,6) * t166 - Ifges(6,6) * t109 - t378;
t364 = t30 / 0.2e1;
t363 = t31 / 0.2e1;
t361 = t377 / 0.2e1;
t359 = t87 / 0.2e1;
t357 = t109 / 0.2e1;
t356 = t110 / 0.2e1;
t186 = t213 * t252;
t187 = t263 * t252;
t128 = -t186 * t253 + t187 * t250;
t355 = t128 / 0.2e1;
t129 = -t186 * t250 - t187 * t253;
t354 = t129 / 0.2e1;
t352 = t277 / 0.2e1;
t351 = -t153 / 0.2e1;
t350 = t153 / 0.2e1;
t349 = t165 / 0.2e1;
t348 = t166 / 0.2e1;
t347 = -t186 / 0.2e1;
t346 = -t187 / 0.2e1;
t345 = -t219 / 0.2e1;
t344 = -t220 / 0.2e1;
t341 = t232 / 0.2e1;
t340 = -t239 / 0.2e1;
t335 = pkin(4) * t220;
t309 = t241 * t251;
t301 = t254 * t223 + t297 * t309;
t67 = -t252 * t291 + t262 * qJD(3) + (-t221 + (qJ(5) * t252 - t209) * t251) * qJD(4) + t301;
t302 = t209 * t292 + t251 * t223;
t307 = t252 * t254;
t75 = (-qJ(5) * qJD(4) - qJD(3) * t241) * t307 + (-qJD(5) * t252 + (-qJ(5) * qJD(3) - qJD(4) * t241) * t255) * t251 + t302;
t34 = t246 * t67 + t248 * t75;
t318 = t227 * mrSges(4,2);
t116 = -t174 * t253 + t175 * t250;
t147 = t213 * t253 - t250 * t263;
t83 = -qJD(6) * t147 - t200 * t253 + t201 * t250;
t313 = t116 - t83;
t117 = -t174 * t250 - t175 * t253;
t146 = -t213 * t250 - t253 * t263;
t82 = qJD(6) * t146 - t200 * t250 - t201 * t253;
t312 = t117 - t82;
t308 = t251 * t252;
t193 = t254 * t209;
t132 = -qJ(5) * t307 + t193 + (-pkin(4) - t309) * t255;
t141 = -qJ(5) * t308 + t155;
t72 = t246 * t132 + t248 * t141;
t288 = mrSges(4,3) * t300;
t305 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t219 + mrSges(5,2) * t220 + t288;
t198 = pkin(4) * t308 + t252 * t241;
t287 = mrSges(4,3) * t299;
t157 = pkin(4) * t379 + t241 * t295;
t243 = -pkin(4) * t254 - pkin(3);
t280 = m(4) * t241 + mrSges(4,3);
t33 = -t246 * t75 + t248 * t67;
t71 = t248 * t132 - t246 * t141;
t275 = mrSges(5,1) * t254 - mrSges(5,2) * t251;
t272 = Ifges(5,1) * t251 + t327;
t270 = Ifges(5,2) * t254 + t328;
t269 = Ifges(5,5) * t251 + Ifges(5,6) * t254;
t58 = -pkin(5) * t255 + t187 * pkin(9) + t71;
t59 = -pkin(9) * t186 + t72;
t21 = -t250 * t59 + t253 * t58;
t22 = t250 * t58 + t253 * t59;
t266 = t120 * t251 - t121 * t254;
t148 = mrSges(5,1) * t279 - mrSges(5,3) * t165;
t149 = -mrSges(5,2) * t279 + mrSges(5,3) * t166;
t265 = -t251 * t148 + t254 * t149;
t169 = -mrSges(5,2) * t239 + mrSges(5,3) * t219;
t170 = mrSges(5,1) * t239 - mrSges(5,3) * t220;
t264 = -t251 * t169 - t254 * t170;
t177 = t189 * qJD(3);
t122 = -t166 * pkin(4) + t177;
t238 = Ifges(5,3) * t279;
t237 = Ifges(6,3) * t279;
t228 = -qJD(3) * mrSges(4,2) + t287;
t179 = pkin(5) * t263 + t243;
t154 = -t255 * t309 + t193;
t145 = pkin(5) * t186 + t198;
t131 = -qJD(3) * t260 - t200 * t252;
t130 = -qJD(3) * t261 + t263 * t293;
t127 = mrSges(6,1) * t239 - mrSges(6,3) * t153;
t126 = -mrSges(6,2) * t239 + mrSges(6,3) * t277;
t115 = pkin(5) * t153 + t335;
t113 = -mrSges(5,1) * t166 + mrSges(5,2) * t165;
t103 = t165 * Ifges(5,1) + t166 * Ifges(5,4) + Ifges(5,5) * t279;
t102 = t165 * Ifges(5,4) + t166 * Ifges(5,2) + Ifges(5,6) * t279;
t101 = -qJD(4) * t155 + t301;
t100 = (-t252 * t296 - t255 * t294) * t241 + t302;
t97 = mrSges(6,1) * t279 - mrSges(6,3) * t110;
t96 = -mrSges(6,2) * t279 + mrSges(6,3) * t109;
t91 = -mrSges(6,1) * t277 + mrSges(6,2) * t153;
t90 = -pkin(5) * t130 + t157;
t78 = Ifges(6,1) * t153 + Ifges(6,5) * t239 + t373;
t74 = mrSges(7,1) * t232 - mrSges(7,3) * t87;
t73 = -mrSges(7,2) * t232 + mrSges(7,3) * t377;
t62 = -t109 * pkin(5) + t122;
t56 = t110 * Ifges(6,1) + t109 * Ifges(6,4) + Ifges(6,5) * t279;
t55 = t110 * Ifges(6,4) + t109 * Ifges(6,2) + Ifges(6,6) * t279;
t49 = -qJD(6) * t129 + t130 * t253 - t131 * t250;
t48 = qJD(6) * t128 + t130 * t250 + t131 * t253;
t41 = -mrSges(7,1) * t377 + mrSges(7,2) * t87;
t26 = -mrSges(7,2) * t279 + mrSges(7,3) * t31;
t25 = mrSges(7,1) * t279 - mrSges(7,3) * t30;
t20 = pkin(9) * t130 + t34;
t19 = pkin(5) * t297 - pkin(9) * t131 + t33;
t7 = t30 * Ifges(7,1) + t31 * Ifges(7,4) + Ifges(7,5) * t279;
t6 = t30 * Ifges(7,4) + t31 * Ifges(7,2) + Ifges(7,6) * t279;
t5 = -qJD(6) * t22 + t19 * t253 - t20 * t250;
t4 = qJD(6) * t21 + t19 * t250 + t20 * t253;
t1 = [m(5) * (t121 * t100 + t120 * t101 + t66 * t154 + t65 * t155) + t130 * t376 + (Ifges(6,5) * t131 + Ifges(6,6) * t130) * t339 + (Ifges(7,1) * t48 + Ifges(7,4) * t49) * t359 + (Ifges(7,4) * t48 + Ifges(7,2) * t49) * t361 + (Ifges(7,4) * t129 + Ifges(7,2) * t128) * t363 + (Ifges(7,1) * t129 + Ifges(7,4) * t128) * t364 + m(7) * (t10 * t4 + t145 * t62 + t2 * t22 + t21 * t3 + t5 * t9 + t89 * t90) + m(6) * (t122 * t198 + t144 * t157 + t15 * t71 + t16 * t72 + t33 * t44 + t34 * t45) + (t10 * t49 + t128 * t2 - t129 * t3 - t48 * t9) * mrSges(7,3) + (t130 * t45 - t131 * t44 + t15 * t187 - t16 * t186) * mrSges(6,3) + (-Ifges(6,1) * t187 - Ifges(6,4) * t186) * t356 + (-Ifges(6,4) * t187 - Ifges(6,2) * t186) * t357 + t122 * (mrSges(6,1) * t186 - mrSges(6,2) * t187) + (Ifges(7,5) * t48 + Ifges(7,6) * t49) * t341 + t56 * t346 + t55 * t347 + (Ifges(6,1) * t131 + Ifges(6,4) * t130) * t350 + (Ifges(6,4) * t131 + Ifges(6,2) * t130) * t352 + t7 * t354 + t6 * t355 + (-t237 / 0.2e1 - t238 / 0.2e1 - t236 / 0.2e1 + t280 * t176 + (0.3e1 / 0.2e1 * t244 + t282 + 0.2e1 * t318 + (-m(4) * t188 + m(5) * t172 + t305) * t241 - t380) * qJD(3) + t365) * t255 + (t102 * t338 + t103 * t337 + t273 * t349 + t271 * t348 + (-t251 * t65 - t254 * t66) * mrSges(5,3) + (mrSges(4,3) + t274) * t177 + (t270 * t345 + t272 * t344 + t269 * t340 + t172 * t275 - t254 * t139 / 0.2e1 + t140 * t338 + t266 * mrSges(5,3)) * qJD(4) + ((t286 * mrSges(4,1) + (t326 / 0.2e1 - t325 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,4)) * t252 + Ifges(7,5) * t354 + Ifges(7,6) * t355 + Ifges(6,5) * t346 + Ifges(6,6) * t347 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - Ifges(7,3) / 0.2e1 - t289) * t255) * qJD(1) - t280 * t189 + t281 - t366) * qJD(3) + (t113 + (m(4) + m(5)) * t177 - t228 * qJD(3)) * t241) * t252 + t49 * t391 + t21 * t25 + t22 * t26 + t48 * t40 / 0.2e1 + t4 * t73 + t5 * t74 + t89 * (-mrSges(7,1) * t49 + mrSges(7,2) * t48) + t90 * t41 + t72 * t96 + t71 * t97 + t34 * t126 + t33 * t127 + t62 * (-mrSges(7,1) * t128 + mrSges(7,2) * t129) + t131 * t78 / 0.2e1 + t144 * (-mrSges(6,1) * t130 + mrSges(6,2) * t131) + t145 * t8 + t154 * t148 + t155 * t149 + t157 * t91 + t100 * t169 + t101 * t170 + t198 * t57; t131 * t126 + t130 * t127 + t128 * t25 + t129 * t26 - t186 * t97 - t187 * t96 + t48 * t73 + t49 * t74 + (-t113 + t375) * t255 + (t264 * qJD(4) + t265) * t252 + ((t169 * t254 - t170 * t251 + t228 - t287) * t255 + (t41 + t91 - t288 + t305) * t252) * qJD(3) + m(7) * (t10 * t48 + t3 * t128 + t2 * t129 - t255 * t62 + t297 * t89 + t9 * t49) + m(6) * (-t122 * t255 + t44 * t130 + t45 * t131 + t144 * t297 - t15 * t186 - t16 * t187) + m(4) * (t176 * t252 - t177 * t255 + (-t188 * t252 + t189 * t255) * qJD(3)) + m(5) * ((-qJD(3) * t266 - t177) * t255 + (qJD(3) * t172 - qJD(4) * t267 + t268) * t252); t381 * t41 + ((m(6) * t144 + t91) * t251 * pkin(4) + (-m(5) * t267 + t264) * pkin(8) + t256) * qJD(4) + t268 * mrSges(5,3) + t102 * t337 + (t82 / 0.2e1 - t117 / 0.2e1) * t40 - t305 * t189 - m(6) * (t144 * t160 + t44 * t60 + t45 * t61) - m(5) * (t120 * t142 + t121 * t143 + t172 * t189) + (t174 / 0.2e1 - t200 / 0.2e1) * t77 + (-Ifges(6,5) * t175 - Ifges(6,6) * t174) * t340 + (-Ifges(6,1) * t175 - Ifges(6,4) * t174) * t351 + (-Ifges(6,4) * t175 - Ifges(6,2) * t174) * t353 + (t175 / 0.2e1 - t201 / 0.2e1) * t78 + t383 * t73 + (t10 * t383 + t179 * t62 + t2 * t69 + t3 * t68 + t381 * t89 + t384 * t9) * m(7) + t384 * t74 + t265 * pkin(8) + m(6) * (t122 * t243 + t133 * t44 + t134 * t45 + t15 * t158 + t159 * t16) + (Ifges(7,1) * t82 + Ifges(7,4) * t83) * t359 + (Ifges(7,1) * t117 + Ifges(7,4) * t116) * t360 + (Ifges(7,4) * t82 + Ifges(7,2) * t83) * t361 + (Ifges(7,4) * t117 + Ifges(7,2) * t116) * t362 + (Ifges(7,4) * t147 + Ifges(7,2) * t146) * t363 + (Ifges(7,1) * t147 + Ifges(7,4) * t146) * t364 + (-mrSges(4,1) - t275) * t177 + (mrSges(7,1) * t313 - mrSges(7,2) * t312) * t89 + (-t10 * t313 + t146 * t2 - t147 * t3 + t312 * t9) * mrSges(7,3) + ((t282 - t318 + t390 + t380) * t255 + ((t329 / 0.2e1 + (t392 + Ifges(4,2) / 0.2e1) * t255) * qJD(1) + t281 + t189 * mrSges(4,3) + t366) * t252 + (Ifges(6,5) * t213 + Ifges(7,5) * t147 - Ifges(6,6) * t263 + Ifges(7,6) * t146 + t269) * t297 / 0.2e1) * qJD(1) + (-mrSges(6,1) * t399 - mrSges(6,2) * t303) * t144 + (-t15 * t213 - t16 * t263 + t303 * t44 + t399 * t45) * mrSges(6,3) + (-Ifges(6,5) * t201 - Ifges(6,6) * t200) * t339 + (-Ifges(6,1) * t201 - Ifges(6,4) * t200) * t350 + (-Ifges(6,4) * t201 - Ifges(6,2) * t200) * t352 + (Ifges(6,1) * t213 - Ifges(6,4) * t263) * t356 + (Ifges(6,4) * t213 - Ifges(6,2) * t263) * t357 + t122 * (mrSges(6,1) * t263 + mrSges(6,2) * t213) - t263 * t55 / 0.2e1 + (Ifges(7,5) * t82 + Ifges(7,6) * t83) * t341 + (Ifges(7,5) * t117 + Ifges(7,6) * t116) * t342 + t270 * t348 + t272 * t349 + (t83 / 0.2e1 - t116 / 0.2e1) * t39 + t251 * t103 / 0.2e1 + m(5) * (-pkin(3) * t177 + pkin(8) * t268) + t68 * t25 + t69 * t26 - pkin(3) * t113 + t397 * t126 + t398 * t127 + t146 * t6 / 0.2e1 + t147 * t7 / 0.2e1 + t62 * (-mrSges(7,1) * t146 + mrSges(7,2) * t147) + t158 * t97 + t159 * t96 - t160 * t91 - t143 * t169 - t142 * t170 - t176 * mrSges(4,2) + t179 * t8 + t213 * t56 / 0.2e1 - t188 * t228 + t243 * t57; (-t144 * t335 - t44 * t50 - t45 * t51 + (t15 * t248 + t16 * t246) * pkin(4)) * m(6) + t396 + (t120 * t219 + t121 * t220) * mrSges(5,3) + (Ifges(6,1) * t277 - t382) * t351 + t237 + t238 + t87 * t391 + t153 * t376 + (-Ifges(6,2) * t153 + t373 + t78) * t353 + (t153 * t45 + t277 * t44) * mrSges(6,3) + (Ifges(5,5) * t219 + Ifges(6,5) * t277 - Ifges(5,6) * t220 - Ifges(6,6) * t153) * t340 - t144 * (mrSges(6,1) * t153 + mrSges(6,2) * t277) + t139 * t343 + (Ifges(5,1) * t219 - t320) * t344 + (-t220 * t91 + t246 * t96 + t248 * t97) * pkin(4) + (-Ifges(5,2) * t220 + t140 + t210) * t345 - t365 + t371 * t73 + t372 * t74 + (t10 * t371 - t115 * t89 + t194 * t3 + t195 * t2 + t372 * t9) * m(7) - t115 * t41 - t51 * t126 - t50 * t127 - t120 * t169 + t121 * t170 + t194 * t25 + t195 * t26 - t172 * (mrSges(5,1) * t220 + mrSges(5,2) * t219); -t277 * t126 + t153 * t127 - t377 * t73 + t87 * t74 + (-t10 * t377 + t87 * t9 + t62) * m(7) + (t153 * t44 - t277 * t45 + t122) * m(6) - t375; t10 * t74 + t39 * t359 - t9 * t73 + t378 + t396;];
tauc  = t1(:);
