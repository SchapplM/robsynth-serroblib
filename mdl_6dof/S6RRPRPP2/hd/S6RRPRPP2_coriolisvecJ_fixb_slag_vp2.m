% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRPP2
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:50:00
% EndTime: 2019-03-09 09:50:26
% DurationCPUTime: 16.66s
% Computational Cost: add. (5961->592), mult. (15555->743), div. (0->0), fcn. (10621->6), ass. (0->260)
t403 = Ifges(7,4) + Ifges(6,5);
t405 = -Ifges(5,4) + t403;
t391 = Ifges(6,4) + Ifges(5,5);
t359 = Ifges(7,5) - t391;
t194 = sin(pkin(9));
t196 = sin(qJ(2));
t299 = cos(pkin(9));
t330 = cos(qJ(2));
t233 = t299 * t330;
t201 = -t194 * t196 + t233;
t163 = t201 * qJD(2);
t150 = qJD(1) * t163;
t197 = cos(qJ(4));
t243 = t299 * t196;
t254 = qJD(1) * t330;
t161 = -qJD(1) * t243 - t194 * t254;
t195 = sin(qJ(4));
t215 = t197 * qJD(2) + t161 * t195;
t396 = qJD(4) * t215;
t85 = t150 * t197 + t396;
t351 = t85 / 0.2e1;
t360 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t404 = t360 * t351;
t132 = qJD(2) * t195 - t161 * t197;
t86 = qJD(4) * t132 + t150 * t195;
t349 = t86 / 0.2e1;
t176 = t194 * t330 + t243;
t162 = t176 * qJD(2);
t149 = qJD(1) * t162;
t343 = -t149 / 0.2e1;
t342 = t149 / 0.2e1;
t389 = Ifges(7,2) + Ifges(6,3);
t388 = Ifges(5,6) - Ifges(6,6);
t387 = Ifges(7,6) - Ifges(6,6);
t386 = -Ifges(5,3) - Ifges(6,2);
t285 = qJD(1) * t196;
t160 = qJD(1) * t233 - t194 * t285;
t156 = qJD(4) - t160;
t183 = (-qJ(3) - pkin(7)) * t196;
t178 = qJD(1) * t183;
t271 = t330 * pkin(7);
t184 = qJ(3) * t330 + t271;
t179 = t184 * qJD(1);
t244 = t299 * t179;
t121 = t178 * t194 + t244;
t401 = qJD(5) * t195 + t121;
t400 = t403 * t132;
t399 = t403 * t215;
t398 = t403 * t197;
t397 = t403 * t195;
t348 = pkin(4) + pkin(5);
t170 = qJD(2) * pkin(2) + t178;
t120 = t194 * t170 + t244;
t114 = qJD(2) * pkin(8) + t120;
t280 = qJD(4) * t197;
t281 = qJD(4) * t195;
t158 = qJD(2) * t183 + qJD(3) * t330;
t139 = t158 * qJD(1);
t159 = -qJD(2) * t184 - t196 * qJD(3);
t199 = qJD(1) * t159;
t84 = t139 * t299 + t194 * t199;
t277 = qJD(1) * qJD(2);
t252 = t196 * t277;
t240 = pkin(2) * t252;
t89 = pkin(3) * t149 - pkin(8) * t150 + t240;
t191 = -pkin(2) * t330 - pkin(1);
t286 = qJD(1) * t191;
t180 = qJD(3) + t286;
t96 = -t160 * pkin(3) + t161 * pkin(8) + t180;
t9 = -t114 * t280 - t195 * t84 + t197 * t89 - t96 * t281;
t1 = -qJ(6) * t85 - qJD(6) * t132 - t149 * t348 - t9;
t395 = -mrSges(7,3) * t1 - t342 * t359 + t349 * t405 + t404;
t8 = -t114 * t281 + t195 * t89 + t197 * t84 + t96 * t280;
t5 = t149 * qJ(5) + t156 * qJD(5) + t8;
t2 = qJ(6) * t86 - qJD(6) * t215 + t5;
t394 = t2 * mrSges(7,3) - t85 * Ifges(5,4) / 0.2e1 + t403 * t351 + (Ifges(5,2) + t389) * t349 + (Ifges(5,6) + t387) * t343;
t393 = -t162 / 0.2e1;
t392 = t163 / 0.2e1;
t384 = -t156 * t387 - t215 * t389 + t400;
t383 = t132 * t391 - t156 * t386 + t215 * t388;
t60 = mrSges(5,1) * t149 - mrSges(5,3) * t85;
t61 = -t149 * mrSges(6,1) + t85 * mrSges(6,2);
t382 = t61 - t60;
t90 = mrSges(6,2) * t215 + mrSges(6,3) * t156;
t316 = mrSges(5,3) * t215;
t92 = -mrSges(5,2) * t156 + t316;
t320 = t90 + t92;
t315 = mrSges(5,3) * t132;
t94 = mrSges(5,1) * t156 - t315;
t95 = -mrSges(6,1) * t156 + mrSges(6,2) * t132;
t319 = t94 - t95;
t301 = t161 * Ifges(4,4);
t373 = Ifges(4,6) * qJD(2);
t380 = t132 * Ifges(7,5) + t160 * Ifges(4,2) - Ifges(7,6) * t215 - t156 * Ifges(7,3) - t301 + t373;
t270 = pkin(2) * t285;
t106 = -pkin(3) * t161 - pkin(8) * t160 + t270;
t164 = t194 * t179;
t122 = t178 * t299 - t164;
t116 = t195 * t122;
t329 = pkin(2) * t194;
t189 = pkin(8) + t329;
t287 = qJ(6) - t189;
t174 = t287 * t197;
t379 = -t116 - (-qJ(6) * t160 - t106) * t197 - t348 * t161 - qJD(4) * t174 - qJD(6) * t195;
t296 = t160 * t195;
t46 = t195 * t106 + t197 * t122;
t37 = -t161 * qJ(5) + t46;
t378 = -qJ(6) * t296 - qJD(6) * t197 + t281 * t287 - t37;
t297 = qJ(5) * t197;
t213 = -t195 * t348 + t297;
t377 = t156 * t213 + t401;
t220 = pkin(4) * t195 - t297;
t376 = t156 * t220 - t401;
t302 = t161 * mrSges(4,3);
t375 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t215 + mrSges(5,2) * t132 - t302;
t374 = Ifges(4,5) * qJD(2);
t372 = t197 * t348;
t40 = -t195 * t114 + t197 * t96;
t28 = qJ(6) * t132 + t40;
t371 = qJD(5) - t28;
t370 = qJD(5) - t40;
t369 = t162 * qJ(5) - qJD(5) * t201;
t368 = -t195 * t388 + t197 * t391;
t367 = t195 * t389 + t398;
t366 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t285) * t271 + t196 * pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t254);
t261 = t176 * t280;
t288 = t195 * t163;
t210 = t261 + t288;
t365 = -t149 * t386 - t388 * t86 + t391 * t85;
t295 = t160 * t197;
t364 = t280 - t295;
t363 = t281 - t296;
t362 = -t195 * t9 + t197 * t8;
t7 = -pkin(4) * t149 - t9;
t361 = t195 * t7 + t197 * t5;
t130 = Ifges(5,4) * t215;
t357 = t132 * t360 - t156 * t359 + t130 - t399;
t230 = -mrSges(7,1) * t195 + mrSges(7,2) * t197;
t231 = mrSges(6,1) * t195 - mrSges(6,3) * t197;
t232 = mrSges(5,1) * t195 + mrSges(5,2) * t197;
t119 = t299 * t170 - t164;
t234 = qJD(2) * pkin(3) + t119;
t205 = qJ(5) * t132 + t234;
t27 = t215 * t348 + qJD(6) + t205;
t42 = -pkin(4) * t215 - t205;
t356 = t27 * t230 + t42 * t231 - t232 * t234;
t312 = Ifges(5,4) * t195;
t355 = t197 * t360 - t312 + t397;
t354 = t9 * mrSges(5,1) - t7 * mrSges(6,1) - t1 * mrSges(7,1) - t8 * mrSges(5,2) + t2 * mrSges(7,2) + t5 * mrSges(6,3);
t350 = -t86 / 0.2e1;
t347 = t215 / 0.2e1;
t346 = -t215 / 0.2e1;
t345 = -t132 / 0.2e1;
t344 = t132 / 0.2e1;
t341 = -t156 / 0.2e1;
t340 = t156 / 0.2e1;
t339 = -t160 / 0.2e1;
t338 = -t161 / 0.2e1;
t337 = t161 / 0.2e1;
t326 = t197 * pkin(4);
t62 = mrSges(7,2) * t149 + mrSges(7,3) * t86;
t64 = -mrSges(6,2) * t86 + mrSges(6,3) * t149;
t323 = t62 + t64;
t69 = -mrSges(6,1) * t215 - mrSges(6,3) * t132;
t70 = mrSges(7,1) * t215 + mrSges(7,2) * t132;
t322 = t69 - t70;
t91 = mrSges(7,2) * t156 - mrSges(7,3) * t215;
t321 = t90 + t91;
t318 = mrSges(4,3) * t149;
t317 = mrSges(4,3) * t150;
t314 = Ifges(3,4) * t196;
t313 = Ifges(5,4) * t132;
t311 = Ifges(5,4) * t197;
t126 = -t299 * t183 + t184 * t194;
t83 = t139 * t194 - t299 * t199;
t306 = t126 * t83;
t303 = t160 * mrSges(4,3);
t300 = t176 * t83;
t41 = t197 * t114 + t195 * t96;
t298 = qJ(5) * t215;
t294 = t163 * t197;
t293 = t176 * t195;
t289 = t195 * qJ(5);
t118 = -pkin(3) * t201 - t176 * pkin(8) + t191;
t127 = t194 * t183 + t184 * t299;
t66 = t195 * t118 + t197 * t127;
t284 = qJD(2) * t196;
t278 = qJD(5) * t197;
t269 = pkin(2) * t284;
t103 = t158 * t299 + t194 * t159;
t107 = pkin(3) * t162 - pkin(8) * t163 + t269;
t268 = t197 * t103 + t195 * t107 + t118 * t280;
t267 = t195 * t103 + t118 * t281 + t127 * t280;
t266 = Ifges(3,4) * t330;
t43 = -qJ(5) * t201 + t66;
t265 = t299 * pkin(2);
t264 = t176 * t281;
t34 = -t86 * mrSges(7,1) + t85 * mrSges(7,2);
t253 = qJD(2) * t330;
t247 = -t281 / 0.2e1;
t245 = t280 / 0.2e1;
t59 = -t149 * mrSges(7,1) - t85 * mrSges(7,3);
t242 = t149 * mrSges(4,1) + t150 * mrSges(4,2);
t45 = t106 * t197 - t116;
t124 = t195 * t127;
t65 = t118 * t197 - t124;
t102 = t158 * t194 - t299 * t159;
t237 = qJD(1) * t253;
t190 = -t265 - pkin(3);
t29 = -qJ(6) * t215 + t41;
t226 = -Ifges(5,2) * t195 + t311;
t221 = Ifges(7,5) * t197 + Ifges(7,6) * t195;
t216 = -qJ(6) * t163 - qJD(6) * t176;
t15 = t107 * t197 - t267;
t214 = Ifges(3,5) * t330 - Ifges(3,6) * t196;
t209 = t264 - t294;
t14 = -t127 * t281 + t268;
t207 = t190 - t289;
t206 = Ifges(7,5) * t85 + Ifges(7,6) * t86 - Ifges(7,3) * t149;
t204 = pkin(1) * (mrSges(3,1) * t196 + mrSges(3,2) * t330);
t203 = t196 * (Ifges(3,1) * t330 - t314);
t202 = (Ifges(3,2) * t330 + t314) * qJD(1);
t200 = qJ(5) * t85 + qJD(5) * t132 - t83;
t192 = Ifges(3,4) * t254;
t173 = t287 * t195;
t172 = t207 - t326;
t168 = Ifges(3,1) * t285 + Ifges(3,5) * qJD(2) + t192;
t167 = Ifges(3,6) * qJD(2) + t202;
t155 = Ifges(4,4) * t160;
t148 = -t207 + t372;
t147 = t156 * qJ(5);
t137 = -qJD(2) * mrSges(4,2) + t303;
t115 = -mrSges(4,1) * t160 - mrSges(4,2) * t161;
t111 = -t161 * Ifges(4,1) + t155 + t374;
t93 = -mrSges(7,1) * t156 - mrSges(7,3) * t132;
t68 = pkin(4) * t132 - t298;
t67 = t176 * t220 + t126;
t63 = -mrSges(5,2) * t149 - mrSges(5,3) * t86;
t55 = Ifges(5,2) * t215 + t156 * Ifges(5,6) + t313;
t48 = t176 * t213 - t126;
t47 = -t132 * t348 + t298;
t44 = pkin(4) * t201 - t65;
t38 = pkin(4) * t161 - t45;
t36 = qJ(6) * t293 + t43;
t35 = mrSges(5,1) * t86 + mrSges(5,2) * t85;
t33 = mrSges(6,1) * t86 - mrSges(6,3) * t85;
t32 = t147 + t41;
t31 = -pkin(4) * t156 + t370;
t30 = t124 + (-qJ(6) * t176 - t118) * t197 + t348 * t201;
t18 = t147 + t29;
t17 = t220 * t163 + (-t278 + (t289 + t326) * qJD(4)) * t176 + t102;
t16 = -t156 * t348 + t371;
t13 = t213 * t163 + (t278 + (-t289 - t372) * qJD(4)) * t176 - t102;
t12 = pkin(4) * t86 - t200;
t11 = -pkin(4) * t162 - t15;
t10 = t14 + t369;
t6 = -t348 * t86 + t200;
t4 = qJ(6) * t261 + (-qJD(4) * t127 - t216) * t195 + t268 + t369;
t3 = qJ(6) * t264 - t348 * t162 + (-t107 + t216) * t197 + t267;
t19 = [(Ifges(4,5) * t392 + Ifges(4,6) * t393 + t214 * qJD(2) / 0.2e1 - t366) * qJD(2) + (t317 + t35) * t126 + (-t162 * t387 - t209 * t403 + t210 * t389) * t346 - t234 * (mrSges(5,1) * t210 - mrSges(5,2) * t209) + (-m(4) * t119 - m(5) * t234 + t375) * t102 - t210 * t55 / 0.2e1 + t111 * t392 + t380 * t393 + (-Ifges(5,4) * t209 - Ifges(5,2) * t210 + Ifges(5,6) * t162) * t347 + (Ifges(4,1) * t163 - Ifges(4,4) * t162) * t338 + (-Ifges(7,5) * t209 + Ifges(7,6) * t210 - Ifges(7,3) * t162) * t341 + t383 * t162 / 0.2e1 + (mrSges(4,2) * t240 + Ifges(4,1) * t150 + t12 * t231 + t226 * t350 + t6 * t230 + t245 * t384 + t247 * t357 + t368 * t342 + t367 * t349 + t355 * t351 + (0.2e1 * Ifges(4,4) + t221) * t343 + (mrSges(6,2) * t7 - mrSges(5,3) * t9 + t395) * t197) * t176 + (t203 - 0.2e1 * t204) * t277 + t191 * t242 + (-t119 * t163 - t120 * t162 + t300) * mrSges(4,3) + (t168 + qJD(1) * (Ifges(3,1) * t196 + t266)) * t253 / 0.2e1 - (t202 + t167) * t284 / 0.2e1 + t357 * t294 / 0.2e1 + t41 * (-mrSges(5,2) * t162 - mrSges(5,3) * t210) + t32 * (-mrSges(6,2) * t210 + mrSges(6,3) * t162) + t27 * (-mrSges(7,1) * t210 - mrSges(7,2) * t209) + m(5) * (t14 * t41 + t15 * t40 + t65 * t9 + t66 * t8 + t306) + m(4) * (t103 * t120 + t306 + t127 * t84 + (t180 + t286) * t269) + m(7) * (t1 * t30 + t13 * t27 + t16 * t3 + t18 * t4 + t2 * t36 + t48 * t6) + m(6) * (t10 * t32 + t11 * t31 + t12 * t67 + t17 * t42 + t43 * t5 + t44 * t7) + (-t5 * mrSges(6,2) - t8 * mrSges(5,3) + t394) * t293 + t16 * (-mrSges(7,1) * t162 + mrSges(7,3) * t209) + t40 * (mrSges(5,1) * t162 + mrSges(5,3) * t209) + t31 * (-mrSges(6,1) * t162 - mrSges(6,2) * t209) + t18 * (mrSges(7,2) * t162 + mrSges(7,3) * t210) + t42 * (mrSges(6,1) * t210 + mrSges(6,3) * t209) + t384 * t288 / 0.2e1 + (-t162 * t386 - t209 * t391 - t210 * t388) * t340 + t232 * t300 + t115 * t269 + (-t359 * t162 - t360 * t209 + t210 * t405) * t344 + (-Ifges(3,2) * t196 + t266) * t237 + (-t365 / 0.2e1 + t206 / 0.2e1 + Ifges(4,4) * t150 - mrSges(4,1) * t240 - Ifges(5,6) * t350 + t84 * mrSges(4,3) + t387 * t349 + t386 * t342 + t359 * t351 - t354 + (0.2e1 * Ifges(4,2) + Ifges(7,3)) * t343) * t201 - t127 * t318 + t48 * t34 + t30 * t59 + t44 * t61 + t36 * t62 + t43 * t64 + t65 * t60 + t66 * t63 + t67 * t33 + t17 * t69 + t13 * t70 + t10 * t90 + t4 * t91 + t14 * t92 + t3 * t93 + t15 * t94 + t11 * t95 + t103 * t137 + t160 * (Ifges(4,4) * t163 - Ifges(4,2) * t162) / 0.2e1 + t180 * (mrSges(4,1) * t162 + mrSges(4,2) * t163); t378 * t91 + (-t1 * t173 + t148 * t6 + t16 * t379 - t174 * t2 + t18 * t378 + t27 * t377) * m(7) + t379 * t93 + t380 * t338 + (t245 - t295 / 0.2e1) * t357 + Ifges(3,5) * t237 + (t311 - t398) * t351 + (t121 * t234 + t190 * t83 - t40 * t45 - t41 * t46) * m(5) + (t366 + (-t203 / 0.2e1 + t204) * qJD(1)) * qJD(1) + (t31 * t364 - t32 * t363 + t361) * mrSges(6,2) + (-t363 * t41 - t364 * t40 + t362) * mrSges(5,3) + (-mrSges(3,1) * t237 + mrSges(3,2) * t252) * pkin(7) + (t368 / 0.2e1 - t221 / 0.2e1) * qJD(4) * t156 + (t301 + t383) * t337 + ((t194 * t84 - t299 * t83) * pkin(2) + t119 * t121 - t120 * t122 - t180 * t270) * m(4) + t312 * t350 + (t247 + t296 / 0.2e1) * t55 - (-Ifges(3,2) * t285 + t168 + t192) * t254 / 0.2e1 + (mrSges(5,2) * t83 + t6 * mrSges(7,2) - t12 * mrSges(6,3) + Ifges(7,5) * t343 + t391 * t342 + t395 + t404) * t195 + (-t83 * mrSges(5,1) - t12 * mrSges(6,1) + t6 * mrSges(7,1) + Ifges(5,2) * t350 - Ifges(7,6) * t343 + t388 * t342 - t389 * t349 - t394) * t197 + t397 * t349 + (t281 / 0.2e1 - t296 / 0.2e1) * t384 + (-Ifges(5,6) * t346 + Ifges(4,2) * t339 + Ifges(7,3) * t340 - t16 * mrSges(7,1) + t40 * mrSges(5,1) - t31 * mrSges(6,1) - t41 * mrSges(5,2) + t32 * mrSges(6,3) + t18 * mrSges(7,2) - t373 / 0.2e1 + t180 * mrSges(4,1) + t387 * t347 + t386 * t341 + t359 * t345) * t161 + (t226 * t346 + Ifges(4,1) * t337 + t221 * t340 - t374 / 0.2e1 - t180 * mrSges(4,2) + t367 * t347 + t368 * t341 + t355 * t345 - t356) * t160 - t375 * t121 + t376 * t69 + t377 * t70 + t119 * t303 + (-t367 / 0.2e1 + t226 / 0.2e1) * t396 + (t355 * t344 + t356) * qJD(4) - Ifges(3,6) * t252 - t115 * t270 - t214 * t277 / 0.2e1 + t167 * t285 / 0.2e1 + (-t16 * t364 + t18 * t363) * mrSges(7,3) - t120 * t302 - t265 * t317 - t318 * t329 - t83 * mrSges(4,1) - t84 * mrSges(4,2) - t37 * t90 - t46 * t92 - t45 * t94 - t38 * t95 + (t111 + t155) * t339 - t122 * t137 + t148 * t34 - Ifges(4,6) * t149 + Ifges(4,5) * t150 + t172 * t33 - t173 * t59 - t174 * t62 + t190 * t35 + (t12 * t172 - t31 * t38 - t32 * t37 + t376 * t42) * m(6) + (-t319 * t280 - t320 * t281 + t382 * t195 + ((-t41 * t195 - t40 * t197) * qJD(4) + t362) * m(5) + (t64 + t63) * t197 + ((-t32 * t195 + t31 * t197) * qJD(4) + t361) * m(6)) * t189; -t160 * t137 + (t322 + t375) * t161 + (-t59 + t156 * (t91 + t320) - t382) * t197 + (t323 + t63 + t156 * (t93 - t319)) * t195 + t242 + (-t1 * t197 - t161 * t27 + t195 * t2 + t156 * (t16 * t195 + t18 * t197)) * m(7) + (t161 * t42 + t195 * t5 - t197 * t7 + t156 * (t195 * t31 + t197 * t32)) * m(6) + (-t234 * t161 + t195 * t8 + t197 * t9 + t156 * (-t195 * t40 + t197 * t41)) * m(5) + (-t119 * t161 - t120 * t160 + t240) * m(4); t365 + (-pkin(4) * t7 + qJ(5) * t5 - t31 * t41 + t32 * t370 - t42 * t68) * m(6) + t55 * t344 + (Ifges(7,5) * t215 + Ifges(7,6) * t132) * t340 + (t132 * t32 - t215 * t31) * mrSges(6,2) + (-t132 * t18 + t16 * t215) * mrSges(7,3) - t27 * (-mrSges(7,1) * t132 + mrSges(7,2) * t215) - t42 * (mrSges(6,1) * t132 - mrSges(6,3) * t215) + t234 * (mrSges(5,1) * t132 + mrSges(5,2) * t215) + (-t132 * t388 + t215 * t391) * t341 + (t2 * qJ(5) - t1 * t348 - t16 * t29 + t18 * t371 - t27 * t47) * m(7) - t348 * t59 + (-Ifges(5,2) * t132 + t130 + t357) * t346 - t206 + (t215 * t360 - t313 + t384 + t400) * t345 + t354 + (t132 * t389 + t399) * t347 + (t316 - t320) * t40 + (t315 + t319) * t41 + t321 * qJD(5) + t323 * qJ(5) - pkin(4) * t61 - t68 * t69 - t47 * t70 - t28 * t91 - t29 * t93; t322 * t132 - t321 * t156 + t59 + t61 + (-t132 * t27 - t156 * t18 + t1) * m(7) + (t132 * t42 - t156 * t32 + t7) * m(6); t215 * t91 + t132 * t93 + 0.2e1 * (t6 / 0.2e1 + t18 * t347 + t16 * t344) * m(7) + t34;];
tauc  = t19(:);
