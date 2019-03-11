% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:45:20
% EndTime: 2019-03-09 10:45:40
% DurationCPUTime: 10.37s
% Computational Cost: add. (9692->554), mult. (23345->723), div. (0->0), fcn. (15897->8), ass. (0->264)
t355 = qJD(2) - qJD(4);
t330 = t355 / 0.2e1;
t222 = sin(qJ(4));
t225 = cos(qJ(4));
t226 = cos(qJ(2));
t278 = qJD(1) * t226;
t223 = sin(qJ(2));
t279 = qJD(1) * t223;
t148 = -t222 * t279 - t225 * t278;
t149 = -t222 * t278 + t225 * t279;
t219 = sin(pkin(10));
t288 = cos(pkin(10));
t362 = t148 * t219 + t288 * t149;
t340 = -t362 / 0.2e1;
t363 = t288 * t148 - t219 * t149;
t341 = -t363 / 0.2e1;
t404 = Ifges(6,4) * t340 + Ifges(6,2) * t341 + Ifges(6,6) * t330;
t221 = sin(qJ(6));
t224 = cos(qJ(6));
t76 = -t221 * t355 + t224 * t362;
t343 = t76 / 0.2e1;
t329 = -t355 / 0.2e1;
t385 = qJD(6) - t363;
t403 = t385 / 0.2e1;
t328 = -t221 / 0.2e1;
t402 = t224 / 0.2e1;
t286 = qJ(5) * t149;
t209 = pkin(7) * t279;
t175 = pkin(8) * t279 - t209;
t241 = qJD(3) - t175;
t227 = -pkin(2) - pkin(3);
t262 = t227 * qJD(2);
t232 = t262 + t241;
t122 = t225 * t232;
t210 = pkin(7) * t278;
t176 = -pkin(8) * t278 + t210;
t218 = qJD(2) * qJ(3);
t150 = t176 + t218;
t83 = -t150 * t222 + t122;
t68 = t83 - t286;
t65 = -pkin(4) * t355 + t68;
t285 = t148 * qJ(5);
t283 = t225 * t150;
t84 = t222 * t232 + t283;
t69 = t84 + t285;
t66 = t288 * t69;
t36 = t219 * t65 + t66;
t32 = -pkin(9) * t355 + t36;
t151 = -qJD(1) * pkin(1) - pkin(2) * t278 - qJ(3) * t279;
t126 = pkin(3) * t278 - t151;
t90 = -pkin(4) * t148 + qJD(5) + t126;
t39 = -pkin(5) * t363 - pkin(9) * t362 + t90;
t10 = -t221 * t32 + t224 * t39;
t11 = t221 * t39 + t224 * t32;
t75 = -t221 * t362 - t224 * t355;
t398 = Ifges(7,3) * t403 + Ifges(7,5) * t343 + t75 * Ifges(7,6) / 0.2e1 + t404;
t401 = t90 * mrSges(6,1) + t10 * mrSges(7,1) - t11 * mrSges(7,2) + t398;
t294 = t219 * t69;
t35 = t288 * t65 - t294;
t316 = t76 * Ifges(7,4);
t25 = t75 * Ifges(7,2) + Ifges(7,6) * t385 + t316;
t74 = Ifges(7,4) * t75;
t26 = t76 * Ifges(7,1) + Ifges(7,5) * t385 + t74;
t359 = t25 * t328 + t26 * t402;
t250 = Ifges(7,5) * t224 - Ifges(7,6) * t221;
t304 = Ifges(7,4) * t224;
t252 = -Ifges(7,2) * t221 + t304;
t305 = Ifges(7,4) * t221;
t253 = Ifges(7,1) * t224 - t305;
t342 = -t385 / 0.2e1;
t344 = -t76 / 0.2e1;
t345 = -t75 / 0.2e1;
t384 = t250 * t342 + t252 * t345 + t253 * t344;
t397 = Ifges(6,5) * t329 + t362 * Ifges(6,1) / 0.2e1 + t363 * Ifges(6,4) / 0.2e1;
t400 = -t35 * mrSges(6,3) - Ifges(6,1) * t340 - Ifges(6,4) * t341 - Ifges(6,5) * t330 + t359 - t384 + t397;
t135 = Ifges(5,4) * t148;
t136 = Ifges(5,4) * t149;
t191 = Ifges(7,2) * t224 + t305;
t192 = Ifges(7,1) * t221 + t304;
t327 = -t224 / 0.2e1;
t334 = t148 / 0.2e1;
t242 = t222 * t223 + t225 * t226;
t111 = t355 * t242;
t97 = t111 * qJD(1);
t243 = t222 * t226 - t223 * t225;
t110 = t355 * t243;
t98 = t110 * qJD(1);
t60 = -t219 * t98 + t288 * t97;
t34 = -qJD(6) * t76 - t221 * t60;
t33 = qJD(6) * t75 + t224 * t60;
t391 = t33 / 0.2e1;
t254 = mrSges(7,1) * t221 + mrSges(7,2) * t224;
t31 = pkin(5) * t355 - t35;
t395 = t31 * t254;
t59 = t219 * t97 + t288 * t98;
t5 = Ifges(7,4) * t33 + Ifges(7,2) * t34 + Ifges(7,6) * t59;
t339 = pkin(7) - pkin(8);
t194 = t339 * t226;
t178 = qJD(2) * t194;
t153 = t222 * t178;
t276 = qJD(4) * t222;
t364 = t339 * t223;
t54 = qJD(4) * t122 - t150 * t276 + t225 * (-qJD(1) * t364 + qJD(3)) * qJD(2) + qJD(1) * t153;
t117 = t225 * t194 + t222 * t364;
t229 = ((-qJD(4) * t227 - qJD(3)) * t222 + t117 * qJD(1)) * qJD(2);
t235 = t222 * t241;
t55 = (-t235 - t283) * qJD(4) + t229;
t6 = t33 * Ifges(7,1) + t34 * Ifges(7,4) + t59 * Ifges(7,5);
t245 = -t97 * qJ(5) - t149 * qJD(5);
t275 = qJD(4) * t225;
t28 = -qJ(5) * t98 + qJD(5) * t148 + t54;
t8 = t288 * t28 + (-qJD(4) * t235 - t150 * t275 + t229 + t245) * t219;
t86 = t148 * Ifges(5,2) - Ifges(5,6) * t355 + t136;
t87 = t149 * Ifges(5,1) - Ifges(5,5) * t355 + t135;
t399 = -(-mrSges(6,3) * t36 - Ifges(7,5) * t344 - Ifges(7,6) * t345 - Ifges(7,3) * t342 + t401 + t404) * t362 - t126 * (mrSges(5,1) * t149 + mrSges(5,2) * t148) - t54 * mrSges(5,2) - t8 * mrSges(6,2) + Ifges(5,5) * t97 + Ifges(6,5) * t60 - Ifges(5,6) * t98 - t6 * t328 - t5 * t327 + t192 * t391 + t34 * t191 / 0.2e1 + t55 * mrSges(5,1) + (Ifges(5,2) * t149 - t135 - t87) * t334 - (Ifges(5,1) * t148 - t136 - t86) * t149 / 0.2e1 + (t395 - t384) * qJD(6);
t369 = -qJD(2) / 0.2e1;
t396 = -mrSges(3,1) - mrSges(4,1);
t249 = -t10 * t224 - t11 * t221;
t353 = t90 * mrSges(6,2) + t395;
t394 = t249 * mrSges(7,3) + t353;
t180 = -qJ(3) * t222 + t225 * t227;
t131 = t225 * qJD(3) + qJD(4) * t180;
t181 = t225 * qJ(3) + t222 * t227;
t132 = -t222 * qJD(3) - qJD(4) * t181;
t106 = -t175 * t222 + t225 * t176;
t239 = t106 + t285;
t107 = t225 * t175 + t222 * t176;
t77 = t107 + t286;
t390 = (t132 - t239) * t288 + (-t131 + t77) * t219;
t388 = -t106 + t132;
t387 = -t107 + t131;
t186 = t210 + t218;
t188 = mrSges(4,2) * t278 + qJD(2) * mrSges(4,3);
t306 = Ifges(3,4) * t223;
t374 = -Ifges(4,5) * t279 / 0.2e1;
t375 = Ifges(4,3) / 0.2e1;
t383 = -(m(4) * t186 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t278 + t188) * pkin(7) - t186 * mrSges(4,2) - t278 * t375 - t374 - (t226 * Ifges(3,2) + t306) * qJD(1) / 0.2e1 + t151 * mrSges(4,1) + (-Ifges(4,6) + Ifges(3,6)) * t369;
t301 = qJD(2) * pkin(2);
t182 = qJD(3) + t209 - t301;
t303 = Ifges(4,5) * t226;
t373 = -Ifges(3,4) * t278 / 0.2e1;
t376 = -Ifges(3,1) / 0.2e1;
t382 = (m(4) * t182 + (mrSges(4,2) + mrSges(3,3)) * t279 + t396 * qJD(2)) * pkin(7) - t151 * mrSges(4,3) + (t223 * Ifges(4,1) - t303) * qJD(1) / 0.2e1 - t279 * t376 - t373 + t182 * mrSges(4,2) - (Ifges(4,4) + Ifges(3,5)) * t369;
t381 = -pkin(5) * t362 + pkin(9) * t363;
t259 = t223 * t262;
t204 = qJ(3) * t278;
t213 = t223 * qJD(3);
t281 = qJD(1) * t213 + qJD(2) * t204;
t113 = qJD(1) * t259 + t281;
t64 = pkin(4) * t98 + t113;
t16 = pkin(5) * t59 - pkin(9) * t60 + t64;
t1 = qJD(6) * t10 + t16 * t221 + t224 * t8;
t2 = -qJD(6) * t11 + t16 * t224 - t221 * t8;
t380 = -t1 * t224 + t2 * t221;
t372 = t252 / 0.2e1;
t231 = qJD(6) * t249 - t380;
t371 = m(7) * t231;
t309 = mrSges(6,1) * t355 - mrSges(7,1) * t75 + mrSges(7,2) * t76 + mrSges(6,3) * t362;
t255 = mrSges(7,1) * t224 - mrSges(7,2) * t221;
t368 = mrSges(6,1) + t255;
t162 = t219 * t225 + t222 * t288;
t360 = t355 * t162;
t19 = mrSges(7,1) * t59 - mrSges(7,3) * t33;
t20 = -mrSges(7,2) * t59 + mrSges(7,3) * t34;
t247 = -t221 * t19 + t224 * t20;
t358 = t148 * t83 + t149 * t84;
t357 = -t10 * t221 + t11 * t224;
t354 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t33 + Ifges(7,6) * t34;
t346 = t59 / 0.2e1;
t338 = pkin(1) * mrSges(3,1);
t337 = pkin(1) * mrSges(3,2);
t116 = -t194 * t222 + t225 * t364;
t238 = qJ(5) * t243 + t116;
t88 = -qJ(5) * t242 + t117;
t50 = t219 * t88 - t238 * t288;
t7 = t219 * t28 - t288 * (t245 + t55);
t336 = t50 * t7;
t332 = t149 / 0.2e1;
t190 = Ifges(7,5) * t221 + Ifges(7,6) * t224;
t331 = -t190 / 0.2e1;
t326 = pkin(4) * t149;
t325 = pkin(4) * t219;
t323 = t10 * mrSges(7,3);
t321 = t11 * mrSges(7,3);
t236 = -t219 * t222 + t225 * t288;
t320 = t236 * t7;
t319 = t59 * mrSges(6,3);
t318 = t60 * mrSges(6,3);
t308 = mrSges(7,3) * t224;
t293 = t221 * mrSges(7,3);
t174 = -pkin(4) + t180;
t109 = t219 * t174 + t288 * t181;
t277 = qJD(2) * t226;
t280 = qJ(3) * t277 + t213;
t274 = qJD(6) * t221;
t273 = qJD(6) * t224;
t183 = -t226 * pkin(2) - t223 * qJ(3) - pkin(1);
t272 = t223 * t301;
t271 = Ifges(3,5) / 0.2e1 + Ifges(4,4) / 0.2e1;
t270 = 0.3e1 / 0.2e1 * Ifges(4,5) - 0.3e1 / 0.2e1 * Ifges(3,4);
t269 = -Ifges(3,6) / 0.2e1 + Ifges(4,6) / 0.2e1;
t268 = m(4) * pkin(7) + mrSges(4,2);
t261 = t288 * pkin(4);
t260 = t59 * mrSges(6,1) + t60 * mrSges(6,2);
t159 = t226 * pkin(3) - t183;
t257 = -t1 * t221 - t2 * t224;
t251 = Ifges(5,5) * t148 - Ifges(5,6) * t149;
t102 = -t219 * t243 + t242 * t288;
t103 = -t219 * t242 - t243 * t288;
t112 = pkin(4) * t242 + t159;
t47 = pkin(5) * t102 - pkin(9) * t103 + t112;
t51 = t219 * t238 + t288 * t88;
t22 = -t221 * t51 + t224 * t47;
t23 = t221 * t47 + t224 * t51;
t119 = mrSges(5,2) * t355 + mrSges(5,3) * t148;
t120 = -mrSges(5,1) * t355 - t149 * mrSges(5,3);
t244 = t119 * t225 - t120 * t222;
t134 = t227 * t279 + t204;
t48 = -mrSges(7,2) * t385 + mrSges(7,3) * t75;
t49 = mrSges(7,1) * t385 - mrSges(7,3) * t76;
t78 = mrSges(6,2) * t355 + mrSges(6,3) * t363;
t240 = -t221 * t49 + t224 * t48 + t78;
t108 = t174 * t288 - t219 * t181;
t121 = t259 + t280;
t177 = qJD(2) * t364;
t70 = -t225 * t177 - t194 * t276 + t275 * t364 + t153;
t100 = t134 - t326;
t72 = pkin(4) * t110 + t121;
t71 = -qJD(4) * t117 + t177 * t222 + t225 * t178;
t230 = -qJ(5) * t111 + qJD(5) * t243 + t71;
t201 = -t261 - pkin(5);
t179 = (qJD(3) - t209) * qJD(2);
t168 = (-t226 * mrSges(4,1) - mrSges(4,3) * t223) * qJD(1);
t142 = t236 * qJD(4);
t141 = t236 * qJD(2);
t137 = t272 - t280;
t125 = qJD(1) * t272 - t281;
t115 = t141 * t224 + t221 * t279;
t114 = -t141 * t221 + t224 * t279;
t105 = -pkin(9) + t109;
t104 = pkin(5) - t108;
t99 = -mrSges(5,1) * t148 + mrSges(5,2) * t149;
t81 = t131 * t288 + t219 * t132;
t62 = -t219 * t110 + t111 * t288;
t61 = t110 * t288 + t111 * t219;
t58 = Ifges(7,3) * t59;
t56 = -mrSges(6,1) * t363 + mrSges(6,2) * t362;
t46 = t326 - t381;
t44 = t219 * t239 + t288 * t77;
t42 = -qJ(5) * t110 - qJD(5) * t242 + t70;
t40 = t100 + t381;
t38 = t288 * t68 - t294;
t37 = t219 * t68 + t66;
t21 = pkin(5) * t61 - pkin(9) * t62 + t72;
t18 = t219 * t230 + t288 * t42;
t17 = t219 * t42 - t230 * t288;
t15 = t221 * t40 + t224 * t44;
t14 = -t221 * t44 + t224 * t40;
t13 = t221 * t46 + t224 * t38;
t12 = -t221 * t38 + t224 * t46;
t9 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t4 = -qJD(6) * t23 - t18 * t221 + t21 * t224;
t3 = qJD(6) * t22 + t18 * t224 + t21 * t221;
t24 = [(-t125 * mrSges(4,3) + (t269 * qJD(2) + (t183 * mrSges(4,1) + t223 * t270 - 0.2e1 * t338) * qJD(1) + t383) * qJD(2)) * t223 + (-t125 * mrSges(4,1) + t268 * t179 + (t271 * qJD(2) + (-0.2e1 * t337 - t183 * mrSges(4,3) - t270 * t226 + (-0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(3,1) + t268 * pkin(7)) * t223) * qJD(1) + t382) * qJD(2)) * t226 + (t250 * t403 + t253 * t343 + t75 * t372 + t359 + t394 + 0.2e1 * t397) * t62 + (t398 + t401) * t61 + (-t110 * t334 + t242 * t98) * Ifges(5,2) + (-t110 * t84 - t111 * t83 - t116 * t97 - t117 * t98 - t242 * t54 + t243 * t55) * mrSges(5,3) + (-t110 * t332 + t111 * t334 - t242 * t97 + t243 * t98) * Ifges(5,4) + t159 * (mrSges(5,1) * t98 + mrSges(5,2) * t97) + m(5) * (t113 * t159 + t116 * t55 + t117 * t54 + t121 * t126 + t70 * t84 + t71 * t83) + m(6) * (t112 * t64 - t17 * t35 + t18 * t36 + t51 * t8 + t72 * t90 + t336) + m(7) * (t1 * t23 + t10 * t4 + t11 * t3 + t17 * t31 + t2 * t22 + t336) + t309 * t17 + (-t35 * t62 - t36 * t61 + t50 * t60 - t51 * t59) * mrSges(6,3) + t113 * (mrSges(5,1) * t242 - mrSges(5,2) * t243) + (t111 * t332 - t243 * t97) * Ifges(5,1) + (Ifges(5,5) * t111 - Ifges(5,6) * t110) * t329 + t126 * (mrSges(5,1) * t110 + mrSges(5,2) * t111) + (t253 * t391 + t34 * t372 + t250 * t346 + t5 * t328 + t6 * t402 + Ifges(6,1) * t60 - Ifges(6,4) * t59 + t64 * mrSges(6,2) + (mrSges(6,3) + t254) * t7 + t257 * mrSges(7,3) + (-mrSges(7,3) * t357 + t191 * t345 + t192 * t344 + t25 * t327 + t255 * t31 + t26 * t328 + t331 * t385) * qJD(6)) * t103 + m(4) * (t125 * t183 + t137 * t151) + t22 * t19 + t23 * t20 + t3 * t48 + t4 * t49 + t50 * t9 + t72 * t56 + t18 * t78 + (-t8 * mrSges(6,3) + t58 / 0.2e1 - Ifges(6,4) * t60 + t64 * mrSges(6,1) + (Ifges(7,3) / 0.2e1 + Ifges(6,2)) * t59 + t354) * t102 - t110 * t86 / 0.2e1 + t111 * t87 / 0.2e1 + t70 * t119 + t71 * t120 + t121 * t99 + t112 * t260 + t137 * t168; t105 * t371 + (-t100 * t90 - t108 * t7 + t109 * t8 + (-t44 + t81) * t36 + t390 * t35) * m(6) + (-t10 * t14 + t104 * t7 - t11 * t15 - t31 * t390 + t357 * t81) * m(7) - t390 * t309 + t387 * t119 + (-t108 * t60 - t109 * t59) * mrSges(6,3) + ((t373 + (t337 + t303 / 0.2e1) * qJD(1) + (-pkin(2) * mrSges(4,2) + (-m(4) * pkin(2) + t396) * pkin(7) + t271) * qJD(2) - t382) * t226 + (t374 + (t338 + t306 / 0.2e1) * qJD(1) + (t375 + t376 - Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t278 + (pkin(7) * mrSges(3,2) - qJ(3) * mrSges(4,2) + t269) * qJD(2) - t383) * t223) * qJD(1) + (t394 + t400) * t363 + (-t180 * t97 - t181 * t98 - t358) * mrSges(5,3) + t251 * t329 + t380 * mrSges(7,3) + (-t126 * t134 + t180 * t55 + t181 * t54 + t387 * t84 + t388 * t83) * m(5) + t388 * t120 + m(4) * (qJ(3) * t179 + qJD(3) * t186) + (Ifges(6,6) + t331) * t59 + ((-t26 / 0.2e1 + t323 - t105 * t49) * t224 + (t25 / 0.2e1 + t321 - t105 * t48) * t221) * qJD(6) + (-m(4) * t151 - t168) * (pkin(2) * t279 - t204) - t399 + t240 * t81 + t247 * t105 - t15 * t48 - t14 * t49 - t44 * t78 - t100 * t56 + t368 * t7 + t104 * t9 - t134 * t99 + t179 * mrSges(4,3) + qJD(3) * t188; -t114 * t49 - t115 * t48 - t141 * t78 - (t9 + t318) * t236 + t244 * qJD(4) + (-t222 * t98 - t225 * t97) * mrSges(5,3) + t240 * t142 + (-t188 - t244) * qJD(2) + (-t319 + (-t221 * t48 - t224 * t49) * qJD(6) + t247) * t162 + (t268 * t277 + (t168 - t56 - t99) * t223) * qJD(1) - m(4) * (qJD(2) * t186 - t151 * t279) - t360 * t309 + (-t10 * t114 - t11 * t115 + t142 * t357 + t231 * t162 - t31 * t360 - t320) * m(7) + (t162 * t8 - t90 * t279 - t320 + (-t141 + t142) * t36 + t360 * t35) * m(6) + (-t126 * t279 + t222 * t54 + t225 * t55 - t355 * (-t222 * t83 + t225 * t84)) * m(5); (-t326 * t90 + t35 * t37 - t36 * t38 + (t219 * t8 - t288 * t7) * pkin(4)) * m(6) - m(7) * (t10 * t12 + t11 * t13) + (t10 * t308 + t11 * t293 - t353 - t400) * t363 + t190 * t346 + t251 * t330 + t1 * t308 - t56 * t326 - t273 * t323 - t319 * t325 - t274 * t321 - t261 * t318 - t2 * t293 + t399 - t13 * t48 - t12 * t49 - Ifges(6,6) * t59 - t38 * t78 + t358 * mrSges(5,3) + t359 * qJD(6) + (m(7) * t201 - t368) * t7 - t83 * t119 + t84 * t120 + (-m(7) * t31 - t309) * t37 + (-t273 * t49 - t274 * t48 + t247 + t371) * (pkin(9) + t325) + t201 * t9; -t363 * t78 - t309 * t362 + (t385 * t48 + t19) * t224 + (-t385 * t49 + t20) * t221 + t260 + (-t31 * t362 + t357 * t385 - t257) * m(7) + (t35 * t362 - t36 * t363 + t64) * m(6); t58 - t31 * (mrSges(7,1) * t76 + mrSges(7,2) * t75) + (Ifges(7,1) * t75 - t316) * t344 + t25 * t343 + (Ifges(7,5) * t75 - Ifges(7,6) * t76) * t342 - t10 * t48 + t11 * t49 + (t10 * t75 + t11 * t76) * mrSges(7,3) + (-Ifges(7,2) * t76 + t26 + t74) * t345 + t354;];
tauc  = t24(:);
