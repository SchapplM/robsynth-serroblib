% Calculate vector of inverse dynamics joint torques for
% S5RRPRP10
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP10_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:25
% EndTime: 2019-12-31 20:09:54
% DurationCPUTime: 16.99s
% Computational Cost: add. (2889->521), mult. (6151->674), div. (0->0), fcn. (3368->6), ass. (0->237)
t386 = Ifges(5,4) + Ifges(6,4);
t387 = Ifges(5,1) + Ifges(6,1);
t368 = Ifges(5,5) + Ifges(6,5);
t385 = Ifges(5,2) + Ifges(6,2);
t367 = Ifges(5,6) + Ifges(6,6);
t161 = cos(qJ(4));
t390 = t386 * t161;
t158 = sin(qJ(4));
t389 = t386 * t158;
t162 = cos(qJ(2));
t266 = qJD(1) * t162;
t103 = -qJD(2) * t158 - t161 * t266;
t159 = sin(qJ(2));
t255 = qJD(1) * qJD(2);
t113 = -t162 * qJDD(1) + t159 * t255;
t114 = qJDD(1) * t159 + t162 * t255;
t262 = qJD(3) * t159;
t283 = qJDD(1) * pkin(1);
t174 = -qJ(3) * t114 - qJD(1) * t262 - t283;
t328 = pkin(2) + pkin(7);
t25 = t113 * t328 + t174;
t260 = qJD(4) * t161;
t261 = qJD(4) * t158;
t101 = t114 * pkin(6);
t229 = qJDD(3) + t101;
t49 = pkin(3) * t114 - qJDD(2) * t328 + t229;
t146 = t159 * qJ(3);
t228 = -pkin(1) - t146;
t66 = (-t162 * t328 + t228) * qJD(1);
t267 = qJD(1) * t159;
t140 = pkin(6) * t267;
t109 = -pkin(3) * t267 - t140;
t341 = -t109 + qJD(3);
t69 = -qJD(2) * t328 + t341;
t3 = t158 * t49 + t161 * t25 + t69 * t260 - t261 * t66;
t264 = qJD(2) * t161;
t180 = t158 * t266 - t264;
t48 = qJD(4) * t180 - qJDD(2) * t158 + t113 * t161;
t2 = qJ(5) * t48 + qJD(5) * t103 + t3;
t388 = t2 * mrSges(6,2);
t370 = -Ifges(4,4) + Ifges(3,5);
t369 = Ifges(4,5) - Ifges(3,6);
t384 = Ifges(5,3) + Ifges(6,3);
t210 = mrSges(5,1) * t158 + mrSges(5,2) * t161;
t246 = m(6) * pkin(4) + mrSges(6,1);
t289 = t161 * mrSges(6,2);
t383 = -t158 * t246 - t210 - t289;
t354 = t158 * t368 + t161 * t367;
t352 = t161 * t385 + t389;
t349 = t158 * t387 + t390;
t157 = -qJ(5) - pkin(7);
t382 = -m(6) * (-pkin(2) + t157) + mrSges(6,3) + m(5) * t328 + mrSges(5,3);
t294 = Ifges(4,6) * t159;
t191 = -t162 * Ifges(4,3) - t294;
t133 = qJD(4) + t267;
t378 = t386 * t103;
t365 = t368 * t133 - t387 * t180 + t378;
t376 = t386 * t180;
t366 = t385 * t103 + t367 * t133 - t376;
t381 = Ifges(4,5) * qJD(2) + qJD(1) * t191 + t158 * t365 + t161 * t366;
t102 = qJDD(4) + t114;
t47 = qJD(4) * t103 + qJDD(2) * t161 + t113 * t158;
t380 = -t385 * t48 / 0.2e1 - t386 * t47 / 0.2e1 - t367 * t102 / 0.2e1;
t160 = sin(qJ(1));
t379 = g(2) * t160;
t208 = t162 * mrSges(4,2) - t159 * mrSges(4,3);
t213 = mrSges(3,1) * t162 - mrSges(3,2) * t159;
t377 = t208 - t213;
t375 = qJD(3) + t140;
t373 = t102 * t368 + t386 * t48 + t387 * t47;
t312 = g(3) * t162;
t372 = -mrSges(6,1) - mrSges(5,1);
t371 = mrSges(5,2) + mrSges(6,2);
t257 = t161 * qJD(5);
t270 = qJ(5) + t328;
t279 = t158 * t159;
t142 = pkin(6) * t266;
t110 = pkin(3) * t266 + t142;
t141 = pkin(2) * t267;
t284 = qJ(3) * t162;
t190 = pkin(7) * t159 - t284;
t73 = qJD(1) * t190 + t141;
t44 = t161 * t110 - t158 * t73;
t364 = t261 * t270 - t257 - (pkin(4) * t162 - qJ(5) * t279) * qJD(1) - t44;
t116 = t270 * t161;
t45 = t158 * t110 + t161 * t73;
t363 = -qJ(5) * t161 * t267 - qJD(4) * t116 - t158 * qJD(5) - t45;
t244 = mrSges(4,1) * t266;
t122 = -qJD(2) * mrSges(4,3) - t244;
t55 = -mrSges(5,1) * t103 - mrSges(5,2) * t180;
t362 = -t122 + t55;
t136 = pkin(4) * t161 + pkin(3);
t361 = pkin(4) * t260 + t136 * t267 + t375;
t327 = pkin(3) + pkin(6);
t125 = t327 * t159;
t150 = t162 * pkin(2);
t269 = t150 + t146;
t215 = pkin(7) * t162 + t269;
t96 = -pkin(1) - t215;
t53 = t158 * t125 + t161 * t96;
t360 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t266 - t122;
t245 = mrSges(4,1) * t267;
t359 = -mrSges(3,3) * t267 - t245 + (mrSges(3,1) - mrSges(4,2)) * qJD(2);
t358 = t159 * t354 + t162 * t384;
t357 = t159 * t352 + t162 * t367;
t356 = t159 * t349 + t162 * t368;
t293 = Ifges(4,6) * t162;
t355 = t159 * (-Ifges(4,2) * t162 + t294) + t162 * (Ifges(4,3) * t159 - t293);
t353 = -t158 * t367 + t161 * t368;
t351 = -t158 * t385 + t390;
t350 = t159 * t369 + t162 * t370;
t348 = t161 * t387 - t389;
t347 = t102 * t384 + t367 * t48 + t368 * t47;
t100 = t113 * pkin(6);
t346 = -t100 * t162 + t101 * t159;
t67 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t100;
t72 = -qJDD(2) * pkin(2) + t229;
t345 = t159 * t72 - t162 * t67;
t27 = t158 * t69 + t161 * t66;
t4 = -qJD(4) * t27 - t158 * t25 + t161 * t49;
t344 = t158 * t3 + t161 * t4;
t1 = pkin(4) * t102 - qJ(5) * t47 + qJD(5) * t180 + t4;
t343 = t1 * t161 + t158 * t2;
t163 = cos(qJ(1));
t342 = g(1) * t163 + t379;
t256 = m(4) + m(5) + m(6);
t139 = Ifges(3,4) * t266;
t340 = Ifges(3,1) * t267 + Ifges(3,5) * qJD(2) + t103 * t367 + t133 * t384 - t180 * t368 + t139;
t339 = -mrSges(5,1) - t246;
t338 = -mrSges(2,1) + t377;
t288 = t162 * mrSges(4,3);
t337 = -t288 + t383 * t162 + (m(4) * pkin(2) - mrSges(4,2) + t382) * t159;
t305 = mrSges(6,2) * t158;
t209 = mrSges(6,1) * t161 - t305;
t211 = mrSges(5,1) * t161 - mrSges(5,2) * t158;
t156 = qJD(2) * qJ(3);
t83 = t156 + t110;
t51 = -pkin(4) * t103 + qJD(5) + t83;
t334 = t51 * t209 + t83 * t211;
t333 = -m(5) * pkin(3) - m(6) * t136 - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t330 = t47 / 0.2e1;
t329 = t48 / 0.2e1;
t326 = t102 / 0.2e1;
t322 = -t180 / 0.2e1;
t315 = pkin(6) * t159;
t148 = t162 * pkin(6);
t304 = mrSges(5,3) * t103;
t303 = mrSges(5,3) * t180;
t302 = mrSges(6,3) * t103;
t301 = mrSges(6,3) * t180;
t300 = Ifges(3,4) * t159;
t299 = Ifges(3,4) * t162;
t287 = t162 * mrSges(6,3);
t280 = t157 * t162;
t278 = t158 * t160;
t277 = t158 * t162;
t276 = t159 * t161;
t275 = t159 * t163;
t274 = t160 * t161;
t273 = t161 * t162;
t272 = t161 * t163;
t271 = t162 * t163;
t126 = t162 * pkin(3) + t148;
t268 = t163 * pkin(1) + t160 * pkin(6);
t265 = qJD(2) * t159;
t263 = qJD(2) * t162;
t259 = qJD(4) * t162;
t258 = qJD(4) * t328;
t243 = t158 * t259;
t11 = -t48 * mrSges(6,1) + t47 * mrSges(6,2);
t230 = -t260 / 0.2e1;
t135 = pkin(4) * t158 + qJ(3);
t26 = -t158 * t66 + t161 * t69;
t227 = -t255 / 0.2e1;
t76 = t114 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t225 = qJ(5) * t162 - t96;
t224 = pkin(2) * t265 - t262;
t212 = mrSges(3,1) * t159 + mrSges(3,2) * t162;
t203 = t162 * Ifges(3,2) + t300;
t192 = -t159 * Ifges(4,2) - t293;
t188 = t158 * t26 - t161 * t27;
t117 = -qJD(2) * pkin(2) + t375;
t120 = -t142 - t156;
t187 = t117 * t162 + t120 * t159;
t186 = t228 - t150;
t17 = qJ(5) * t180 + t26;
t185 = pkin(1) * t212;
t85 = t159 * t272 - t278;
t87 = t158 * t163 + t159 * t274;
t112 = t327 * t263;
t63 = qJD(2) * t190 + t224;
t13 = t158 * t112 + t125 * t260 + t161 * t63 - t261 * t96;
t84 = t186 * qJD(1);
t184 = t84 * (-mrSges(4,2) * t159 - t288);
t183 = t159 * (Ifges(3,1) * t162 - t300);
t177 = t159 * t264 + t243;
t176 = t158 * t265 - t161 * t259;
t59 = -mrSges(6,2) * t133 + t302;
t60 = -mrSges(5,2) * t133 + t304;
t61 = mrSges(6,1) * t133 + t301;
t62 = mrSges(5,1) * t133 + t303;
t175 = (t59 + t60) * t161 + (-t61 - t62) * t158;
t50 = -pkin(3) * t113 - t67;
t167 = -qJD(4) * t188 + t344;
t118 = -pkin(1) - t269;
t115 = t270 * t158;
t111 = t327 * t265;
t108 = -qJ(3) * t266 + t141;
t107 = t208 * qJD(1);
t106 = t161 * t125;
t95 = t211 * t162;
t92 = t161 * t112;
t88 = -t159 * t278 + t272;
t86 = t158 * t275 + t274;
t82 = Ifges(4,4) * qJD(2) + qJD(1) * t192;
t79 = Ifges(3,6) * qJD(2) + qJD(1) * t203;
t78 = pkin(4) * t273 + t126;
t77 = -qJ(3) * t263 + t224;
t75 = mrSges(4,1) * t113 - qJDD(2) * mrSges(4,3);
t56 = -pkin(4) * t243 + (-pkin(6) - t136) * t265;
t54 = -mrSges(6,1) * t103 - mrSges(6,2) * t180;
t52 = -t158 * t96 + t106;
t46 = pkin(2) * t113 + t174;
t36 = -qJ(5) * t273 + t53;
t29 = pkin(4) * t159 + t158 * t225 + t106;
t23 = -mrSges(5,2) * t102 + mrSges(5,3) * t48;
t22 = -mrSges(6,2) * t102 + mrSges(6,3) * t48;
t21 = mrSges(5,1) * t102 - mrSges(5,3) * t47;
t20 = mrSges(6,1) * t102 - mrSges(6,3) * t47;
t18 = qJ(5) * t103 + t27;
t16 = pkin(4) * t133 + t17;
t15 = -pkin(4) * t48 + qJDD(5) + t50;
t14 = -qJD(4) * t53 - t158 * t63 + t92;
t12 = -mrSges(5,1) * t48 + mrSges(5,2) * t47;
t6 = qJ(5) * t177 - t162 * t257 + t13;
t5 = pkin(4) * t263 + t92 + t225 * t260 + (-qJ(5) * t265 - qJD(4) * t125 + qJD(5) * t162 - t63) * t158;
t7 = [t4 * (mrSges(5,1) * t159 + mrSges(5,3) * t277) + t114 * t299 / 0.2e1 + (-Ifges(3,4) * t113 + Ifges(3,5) * qJDD(2) + t347) * t159 / 0.2e1 + m(5) * (-t111 * t83 + t126 * t50 + t13 * t27 + t14 * t26 + t3 * t53 + t4 * t52) + (t159 * t368 - t162 * t349) * t330 + (t159 * t370 - t162 * t369) * qJDD(2) / 0.2e1 + (-m(3) * t268 + t372 * t86 - t371 * t85 + (-m(5) * pkin(7) - mrSges(5,3)) * t271 - t256 * (pkin(2) * t271 + qJ(3) * t275 + t268) + t333 * t160 + (-m(6) * (pkin(4) * t279 - t280) - t287 + t338) * t163) * g(2) + t345 * mrSges(4,1) + t15 * t209 * t162 - t373 * t277 / 0.2e1 + (t183 + t162 * (-Ifges(3,2) * t159 + t299)) * t255 / 0.2e1 + (t381 / 0.2e1 - t79 / 0.2e1 - t360 * pkin(6) + t120 * mrSges(4,1)) * t265 + t46 * t208 - t113 * t203 / 0.2e1 + t113 * t191 / 0.2e1 + t162 * (Ifges(3,4) * t114 - Ifges(3,2) * t113 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t162 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t114 + Ifges(4,3) * t113) / 0.2e1 - t114 * t192 / 0.2e1 - t159 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t114 + Ifges(4,6) * t113) / 0.2e1 + t118 * (-mrSges(4,2) * t113 - mrSges(4,3) * t114) + t126 * t12 + t213 * t283 + t273 * t380 - t111 * t55 - pkin(1) * (mrSges(3,1) * t113 + mrSges(3,2) * t114) + t77 * t107 + t50 * t95 + t78 * t11 + t6 * t59 + t13 * t60 + t5 * t61 + t14 * t62 + t52 * t21 + t53 * t23 + t56 * t54 + t29 * t20 + t36 * t22 + (t372 * t88 + t371 * t87 + (m(3) * pkin(1) - m(4) * t186 - m(6) * (-t135 * t159 - pkin(1)) - m(5) * t228 + t382 * t162 - t338) * t160 + ((-m(3) - t256) * pkin(6) + t333) * t163) * g(1) + t350 * qJD(2) ^ 2 / 0.2e1 + t355 * t227 + (qJD(2) * t356 - t259 * t348) * t322 + (qJD(2) * t357 - t259 * t351) * t103 / 0.2e1 + (qJD(2) * t358 - t259 * t353) * t133 / 0.2e1 + (t1 * t277 - t16 * t176 + t177 * t18 - t2 * t273) * mrSges(6,3) - t159 * t388 + (t159 * t384 - t354 * t162) * t326 + t365 * t162 * t230 - t185 * t255 + (-qJDD(2) * mrSges(3,2) - t75) * t148 + t114 * t159 * Ifges(3,1) + qJD(2) * t184 + t366 * t243 / 0.2e1 + (t159 * t367 - t162 * t352) * t329 + (-qJDD(2) * mrSges(3,1) + t76) * t315 + t3 * (-mrSges(5,2) * t159 - mrSges(5,3) * t273) - t26 * mrSges(5,3) * t176 + m(6) * (t1 * t29 + t15 * t78 + t16 * t5 + t18 * t6 + t2 * t36 + t51 * t56) + Ifges(2,3) * qJDD(1) + m(4) * (t118 * t46 + t77 * t84 + (qJD(2) * t187 + t345) * pkin(6)) + (-t113 * t148 + t114 * t315 + t346) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t346) + (t340 / 0.2e1 - t82 / 0.2e1 - t18 * mrSges(6,2) - t27 * mrSges(5,2) + t16 * mrSges(6,1) + t26 * mrSges(5,1) - t359 * pkin(6) + t117 * mrSges(4,1)) * t263 + t51 * (-mrSges(6,1) * t177 + mrSges(6,2) * t176) + t83 * (-mrSges(5,1) * t177 + mrSges(5,2) * t176) + t27 * mrSges(5,3) * t177 + t1 * mrSges(6,1) * t159; t369 * t113 + t370 * t114 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2) + (mrSges(6,1) * t15 - t328 * t23 + t258 * t62 + t380) * t158 + t359 * t142 + t360 * t140 - (-Ifges(3,2) * t267 + t139 + t340) * t266 / 0.2e1 + (-qJ(3) * t256 * t271 + t163 * t337) * g(1) + t15 * t289 + t82 * t266 / 0.2e1 + t50 * t210 + t135 * t11 - t115 * t22 - t116 * t20 + t361 * t54 + t362 * qJD(3) + (-t183 / 0.2e1 + t185 + t355 / 0.2e1) * qJD(1) ^ 2 - (t103 * t357 + t133 * t358 - t180 * t356) * qJD(1) / 0.2e1 - (t352 * t103 + t354 * t133 - t180 * t349) * qJD(4) / 0.2e1 + (t50 * qJ(3) - t167 * t328 - t26 * t44 - t27 * t45 + t341 * t83) * m(5) + (-t328 * t21 - t258 * t60 + t373 / 0.2e1) * t161 + (-t256 * t284 + t337) * t379 - t117 * t244 - t109 * t55 - t108 * t107 + t100 * mrSges(3,2) - t101 * mrSges(3,1) - t67 * mrSges(4,3) + t72 * mrSges(4,2) - pkin(2) * t76 - t45 * t60 - t44 * t62 + (-pkin(2) * t72 - qJ(3) * t67 - qJD(3) * t120 - t108 * t84) * m(4) - t381 * t267 / 0.2e1 - t120 * t245 + (-m(4) * t187 * pkin(6) - t184 - t18 * (-mrSges(6,2) * t162 + mrSges(6,3) * t276) - t27 * (-mrSges(5,2) * t162 + mrSges(5,3) * t276) - t16 * (mrSges(6,1) * t162 - mrSges(6,3) * t279) - t26 * (mrSges(5,1) * t162 - mrSges(5,3) * t279)) * qJD(1) + (-t75 + t12) * qJ(3) + (-m(4) * t269 - m(6) * (t269 - t280) - t287 - m(5) * t215 + t383 * t159 + t377) * g(3) + t348 * t330 + t350 * t227 + t351 * t329 + t353 * t326 + t363 * t59 + t364 * t61 + (-t1 * t116 - t115 * t2 + t135 * t15 + t16 * t364 + t18 * t363 + t361 * t51) * m(6) - t365 * t261 / 0.2e1 + t366 * t230 + (t79 / 0.2e1 + t334) * t267 + t334 * qJD(4) + (t26 * t261 - t260 * t27 - t312 - t344) * mrSges(5,3) + t342 * t212 + (t16 * t261 - t18 * t260 - t343) * mrSges(6,3); (t20 + t21) * t161 + (t22 + t23) * t158 + (-t54 - t362) * qJD(2) + t256 * t312 + t175 * qJD(4) + ((t107 + t175) * qJD(1) - t342 * t256) * t159 + t76 + (-qJD(2) * t51 + t343 - t133 * (t158 * t16 - t161 * t18)) * m(6) + (-qJD(2) * t83 - t188 * t267 + t167) * m(5) + (qJD(2) * t120 + t267 * t84 + t72) * m(4); -t83 * (-mrSges(5,1) * t180 + mrSges(5,2) * t103) - t51 * (-mrSges(6,1) * t180 + mrSges(6,2) * t103) - (t103 * t368 + t180 * t367) * t133 / 0.2e1 - (t385 * t180 + t365 + t378) * t103 / 0.2e1 + (t387 * t103 + t376) * t180 / 0.2e1 + (-m(6) * (-t16 + t17) - t301 + t61) * t18 + (-t303 + t62) * t27 - (-t161 * t246 + t305) * t312 + (t339 * t87 - t371 * t88) * g(2) + (t339 * t85 + t371 * t86) * g(1) + t16 * t302 + t246 * t1 + g(3) * t95 - t17 * t59 - t3 * mrSges(5,2) + t4 * mrSges(5,1) - t388 + t366 * t322 + t347 + (t304 - t60) * t26 + (t20 - (-m(6) * t51 - t54) * t180) * pkin(4); -t103 * t59 - t180 * t61 + (-g(3) * t159 - t18 * t103 - t16 * t180 - t342 * t162 + t15) * m(6) + t11;];
tau = t7;
