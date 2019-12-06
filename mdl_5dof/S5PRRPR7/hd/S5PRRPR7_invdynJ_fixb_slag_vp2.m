% Calculate vector of inverse dynamics joint torques for
% S5PRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:02
% EndTime: 2019-12-05 16:35:45
% DurationCPUTime: 15.07s
% Computational Cost: add. (3666->557), mult. (8775->800), div. (0->0), fcn. (6601->12), ass. (0->268)
t197 = sin(qJ(2));
t192 = sin(pkin(5));
t281 = qJD(1) * t192;
t258 = t197 * t281;
t163 = qJD(2) * pkin(7) + t258;
t196 = sin(qJ(3));
t155 = t196 * t163;
t199 = cos(qJ(3));
t194 = cos(pkin(5));
t280 = qJD(1) * t194;
t107 = t199 * t280 - t155;
t379 = qJD(4) - t107;
t221 = pkin(3) * t196 - qJ(4) * t199;
t160 = t221 * qJD(2);
t191 = sin(pkin(10));
t193 = cos(pkin(10));
t375 = -t191 * t160 + t379 * t193;
t200 = cos(qJ(2));
t283 = t199 * t200;
t115 = (t191 * t197 + t193 * t283) * t192;
t141 = qJD(3) * t221 - qJD(4) * t196;
t378 = qJD(1) * t115 - t191 * t141;
t368 = m(6) + m(5);
t277 = qJD(3) * t196;
t377 = -(-pkin(7) * t193 + pkin(8)) * t277 + t378;
t272 = t196 * qJD(2);
t376 = -pkin(8) * t272 + t375;
t232 = pkin(4) * t191 - pkin(8) * t193;
t218 = pkin(7) + t232;
t257 = t200 * t281;
t234 = t196 * t257;
t276 = qJD(3) * t199;
t374 = t218 * t276 - t234;
t270 = qJD(2) * qJD(3);
t162 = qJDD(2) * t196 + t199 * t270;
t120 = qJDD(3) * t193 - t162 * t191;
t373 = -t120 / 0.2e1;
t121 = qJDD(3) * t191 + t162 * t193;
t372 = -t121 / 0.2e1;
t161 = -t199 * qJDD(2) + t196 * t270;
t367 = -t161 / 0.2e1;
t362 = qJD(3) / 0.2e1;
t157 = qJD(3) * t193 - t191 * t272;
t158 = qJD(3) * t191 + t193 * t272;
t265 = mrSges(4,3) * t272;
t353 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t157 + mrSges(5,2) * t158 + t265;
t279 = qJD(2) * t192;
t251 = qJD(1) * t279;
t179 = t200 * t251;
t269 = qJDD(1) * t192;
t137 = t197 * t269 + t179;
t371 = qJDD(2) * pkin(7) + qJD(3) * t280 + t137;
t230 = -mrSges(5,1) * t193 + mrSges(5,2) * t191;
t195 = sin(qJ(5));
t198 = cos(qJ(5));
t228 = mrSges(6,1) * t198 - mrSges(6,2) * t195;
t343 = -m(6) * pkin(4) - t228;
t345 = -m(6) * pkin(8) - mrSges(6,3);
t370 = pkin(3) * t368 - t191 * t345 - t193 * t343 + mrSges(4,1) - t230;
t323 = t161 / 0.2e1;
t325 = t121 / 0.2e1;
t369 = Ifges(5,1) * t325 + Ifges(5,5) * t323;
t271 = t199 * qJD(2);
t99 = -t158 * t195 - t198 * t271;
t28 = qJD(5) * t99 + t121 * t198 + t161 * t195;
t333 = t28 / 0.2e1;
t213 = -t158 * t198 + t195 * t271;
t29 = qJD(5) * t213 - t121 * t195 + t161 * t198;
t332 = t29 / 0.2e1;
t111 = qJDD(5) - t120;
t327 = t111 / 0.2e1;
t326 = t120 / 0.2e1;
t366 = t162 / 0.2e1;
t67 = mrSges(5,1) * t161 - mrSges(5,3) * t121;
t8 = -mrSges(6,1) * t29 + mrSges(6,2) * t28;
t365 = t8 - t67;
t322 = pkin(3) * t199;
t168 = -qJ(4) * t196 - pkin(2) - t322;
t110 = qJD(2) * t168 - t257;
t256 = t196 * t280;
t108 = t163 * t199 + t256;
t96 = qJD(3) * qJ(4) + t108;
t41 = t191 * t110 + t193 * t96;
t32 = -pkin(8) * t271 + t41;
t92 = -qJD(3) * pkin(3) + t379;
t38 = -pkin(4) * t157 - pkin(8) * t158 + t92;
t11 = -t195 * t32 + t198 * t38;
t364 = t11 * mrSges(6,1);
t12 = t195 * t38 + t198 * t32;
t363 = t12 * mrSges(6,2);
t361 = mrSges(3,2) - mrSges(4,3);
t131 = t218 * t196;
t288 = t193 * t199;
t117 = pkin(7) * t288 + t191 * t168;
t98 = -pkin(8) * t199 + t117;
t49 = t131 * t195 + t198 * t98;
t360 = -qJD(5) * t49 + t377 * t195 + t374 * t198;
t48 = t131 * t198 - t195 * t98;
t359 = qJD(5) * t48 + t374 * t195 - t377 * t198;
t165 = -pkin(4) * t193 - pkin(8) * t191 - pkin(3);
t300 = qJ(4) * t193;
t112 = t165 * t198 - t195 * t300;
t65 = t256 + (qJD(2) * t232 + t163) * t199;
t358 = qJD(5) * t112 - t195 * t65 + t376 * t198;
t113 = t165 * t195 + t198 * t300;
t357 = -qJD(5) * t113 - t376 * t195 - t198 * t65;
t231 = mrSges(4,1) * t199 - mrSges(4,2) * t196;
t356 = -mrSges(3,1) - t231;
t267 = pkin(7) * t277;
t355 = -t193 * t267 - t378;
t119 = -mrSges(5,1) * t271 - mrSges(5,3) * t158;
t43 = -mrSges(6,1) * t99 - mrSges(6,2) * t213;
t304 = t119 - t43;
t52 = -t120 * mrSges(5,1) + t121 * mrSges(5,2);
t354 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t162 + t52;
t315 = Ifges(5,4) * t193;
t224 = -Ifges(5,2) * t191 + t315;
t316 = Ifges(5,4) * t191;
t226 = Ifges(5,1) * t193 - t316;
t352 = t157 * (Ifges(5,6) * t196 + t199 * t224) + t158 * (Ifges(5,5) * t196 + t199 * t226);
t268 = qJDD(1) * t194;
t260 = t196 * t268 + t199 * t371;
t36 = -t163 * t277 + t260;
t37 = -t163 * t276 - t196 * t371 + t199 * t268;
t350 = -t196 * t37 + t199 * t36;
t30 = qJDD(3) * qJ(4) + (qJD(4) - t155) * qJD(3) + t260;
t178 = t197 * t251;
t136 = t200 * t269 - t178;
t122 = -qJDD(2) * pkin(2) - t136;
t42 = pkin(3) * t161 - qJ(4) * t162 - qJD(4) * t272 + t122;
t13 = -t191 * t30 + t193 * t42;
t14 = t191 * t42 + t193 * t30;
t349 = -t13 * t191 + t14 * t193;
t348 = -t11 * t195 + t12 * t198;
t154 = qJD(5) - t157;
t318 = Ifges(4,4) * t196;
t225 = Ifges(4,2) * t199 + t318;
t344 = Ifges(4,6) * t362 + qJD(2) * t225 / 0.2e1 - t158 * Ifges(5,5) / 0.2e1 - t157 * Ifges(5,6) / 0.2e1 + Ifges(5,3) * t271 / 0.2e1;
t342 = mrSges(5,2) + t345;
t164 = -qJD(2) * pkin(2) - t257;
t229 = mrSges(5,1) * t191 + mrSges(5,2) * t193;
t292 = t191 * t199;
t40 = t110 * t193 - t191 * t96;
t341 = -t199 * t92 * t229 - t41 * (-mrSges(5,2) * t196 - mrSges(5,3) * t292) - t40 * (mrSges(5,1) * t196 - mrSges(5,3) * t288) - t164 * (mrSges(4,1) * t196 + mrSges(4,2) * t199);
t33 = -qJDD(3) * pkin(3) + qJDD(4) - t37;
t17 = -pkin(4) * t120 - pkin(8) * t121 + t33;
t7 = pkin(8) * t161 + t14;
t1 = qJD(5) * t11 + t17 * t195 + t198 * t7;
t2 = -qJD(5) * t12 + t17 * t198 - t195 * t7;
t340 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t227 = t195 * mrSges(6,1) + t198 * mrSges(6,2);
t339 = -t368 * qJ(4) + mrSges(4,2) - mrSges(5,3) - t227;
t189 = Ifges(4,4) * t271;
t338 = t191 * (-Ifges(6,5) * t213 + t99 * Ifges(6,6) + t154 * Ifges(6,3)) + t193 * (t158 * Ifges(5,1) + t157 * Ifges(5,4) - Ifges(5,5) * t271) + Ifges(4,1) * t272 + Ifges(4,5) * qJD(3) + t189;
t285 = t196 * t198;
t287 = t195 * t196;
t337 = t287 * mrSges(6,1) + t285 * mrSges(6,2) - t356;
t336 = -mrSges(5,1) + t343;
t3 = Ifges(6,5) * t28 + Ifges(6,6) * t29 + Ifges(6,3) * t111;
t335 = Ifges(6,5) * t333 + Ifges(6,6) * t332 + Ifges(6,3) * t327 + t3 / 0.2e1 + Ifges(5,4) * t372 + Ifges(5,2) * t373 + Ifges(5,6) * t367 + t340;
t201 = qJD(2) ^ 2;
t334 = Ifges(6,1) * t333 + Ifges(6,4) * t332 + Ifges(6,5) * t327;
t331 = Ifges(5,4) * t326 + t369;
t330 = -t99 / 0.2e1;
t329 = t213 / 0.2e1;
t328 = -t213 / 0.2e1;
t324 = -t154 / 0.2e1;
t317 = Ifges(4,4) * t199;
t314 = Ifges(6,4) * t213;
t313 = Ifges(6,4) * t195;
t312 = Ifges(6,4) * t198;
t307 = t196 * t33;
t303 = cos(pkin(9));
t302 = sin(pkin(9));
t299 = t141 * t193;
t239 = t302 * t197;
t240 = t303 * t200;
t142 = -t194 * t240 + t239;
t298 = t142 * t196;
t238 = t302 * t200;
t241 = t303 * t197;
t144 = t194 * t238 + t241;
t297 = t144 * t196;
t296 = t168 * t193;
t295 = t191 * t195;
t293 = t191 * t198;
t291 = t192 * t197;
t290 = t192 * t200;
t286 = t195 * t199;
t284 = t198 * t199;
t282 = pkin(2) * t290 + pkin(7) * t291;
t278 = qJD(2) * t197;
t275 = qJD(4) * t191;
t273 = qJD(5) * t191;
t264 = mrSges(4,3) * t271;
t262 = t196 * t290;
t261 = t192 * t283;
t259 = pkin(7) * t191 + pkin(4);
t254 = t192 * t278;
t253 = t200 * t279;
t252 = t191 * t276;
t248 = -t273 / 0.2e1;
t143 = t194 * t241 + t238;
t245 = -t142 * pkin(2) + pkin(7) * t143;
t145 = -t194 * t239 + t240;
t244 = -t144 * pkin(2) + pkin(7) * t145;
t243 = t192 * t303;
t242 = t192 * t302;
t237 = -t270 / 0.2e1;
t235 = pkin(3) * t261 + qJ(4) * t262 + t282;
t233 = t199 * t257;
t223 = Ifges(4,5) * t199 - Ifges(4,6) * t196;
t222 = Ifges(5,5) * t193 - Ifges(5,6) * t191;
t55 = -mrSges(6,2) * t154 + mrSges(6,3) * t99;
t56 = mrSges(6,1) * t154 + mrSges(6,3) * t213;
t219 = -t195 * t56 + t198 * t55;
t147 = -t194 * t199 + t196 * t291;
t148 = t194 * t196 + t199 * t291;
t84 = t148 * t193 - t191 * t290;
t46 = t147 * t198 - t195 * t84;
t47 = t147 * t195 + t198 * t84;
t57 = -t107 * t191 + t160 * t193;
t215 = t193 * t284 + t287;
t151 = t193 * t285 - t286;
t214 = -t193 * t286 + t285;
t150 = -t193 * t287 - t284;
t211 = t196 * (Ifges(4,1) * t199 - t318);
t85 = t143 * t196 + t199 * t243;
t87 = t145 * t196 - t199 * t242;
t210 = -g(1) * t87 - g(2) * t85 - g(3) * t147;
t203 = t199 * (Ifges(5,3) * t196 + t199 * t222);
t172 = -qJD(3) * mrSges(4,2) + t264;
t159 = t231 * qJD(2);
t138 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t161;
t133 = t215 * qJD(2);
t132 = t214 * qJD(2);
t118 = mrSges(5,2) * t271 + mrSges(5,3) * t157;
t116 = -pkin(7) * t292 + t296;
t101 = t191 * t233 - t193 * t258;
t97 = t199 * t259 - t296;
t95 = Ifges(6,4) * t99;
t91 = mrSges(4,1) * t161 + mrSges(4,2) * t162;
t90 = -qJD(3) * t147 + t199 * t253;
t89 = qJD(3) * t148 + t196 * t253;
t88 = t145 * t199 + t196 * t242;
t86 = t143 * t199 - t196 * t243;
t83 = t148 * t191 + t193 * t290;
t81 = t191 * t267 + t299;
t73 = -t259 * t277 - t299;
t72 = qJD(3) * t214 - qJD(5) * t151;
t71 = qJD(3) * t215 + qJD(5) * t150;
t69 = t158 * Ifges(5,4) + t157 * Ifges(5,2) - Ifges(5,6) * t271;
t66 = -mrSges(5,2) * t161 + mrSges(5,3) * t120;
t54 = t191 * t254 + t193 * t90;
t53 = t191 * t90 - t193 * t254;
t50 = -pkin(4) * t272 - t57;
t45 = t144 * t191 + t193 * t88;
t44 = t142 * t191 + t193 * t86;
t31 = pkin(4) * t271 - t40;
t24 = -Ifges(6,1) * t213 + Ifges(6,5) * t154 + t95;
t23 = Ifges(6,2) * t99 + Ifges(6,6) * t154 - t314;
t21 = -mrSges(6,2) * t111 + mrSges(6,3) * t29;
t20 = mrSges(6,1) * t111 - mrSges(6,3) * t28;
t10 = qJD(5) * t46 + t195 * t89 + t198 * t54;
t9 = -qJD(5) * t47 - t195 * t54 + t198 * t89;
t6 = -pkin(4) * t161 - t13;
t4 = t28 * Ifges(6,4) + t29 * Ifges(6,2) + t111 * Ifges(6,6);
t5 = [m(2) * qJDD(1) + t10 * t55 + t54 * t118 + t148 * t138 + t90 * t172 + t46 * t20 + t47 * t21 + t9 * t56 + t84 * t66 + t353 * t89 + t365 * t83 - t304 * t53 + t354 * t147 + (-m(2) - m(3) - m(4) - t368) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t201 - t91) * t200 + (-mrSges(3,1) * t201 - mrSges(3,2) * qJDD(2) - qJD(2) * t159) * t197) * t192 + m(5) * (-t13 * t83 + t14 * t84 + t147 * t33 - t40 * t53 + t41 * t54 + t89 * t92) + m(6) * (t1 * t47 + t10 * t12 + t11 * t9 + t2 * t46 + t31 * t53 + t6 * t83) + m(3) * (qJDD(1) * t194 ^ 2 + (t136 * t200 + t137 * t197) * t192) + m(4) * (-t107 * t89 + t108 * t90 - t147 * t37 + t148 * t36 + (-t122 * t200 + t164 * t278) * t192); t81 * t119 + t116 * t67 + t117 * t66 - pkin(2) * t91 + t97 * t8 + t71 * t24 / 0.2e1 + t31 * (-mrSges(6,1) * t72 + mrSges(6,2) * t71) + t72 * t23 / 0.2e1 + t73 * t43 + t48 * t20 + t49 * t21 + t359 * t55 + (t1 * t49 + t2 * t48 + t6 * t97 + (-t101 + t73) * t31 + t359 * t12 + t360 * t11) * m(6) + t360 * t56 + (-m(4) * t282 - m(6) * (pkin(4) * t115 + t235) - (t115 * t198 + t195 * t262) * mrSges(6,1) - (-t115 * t195 + t198 * t262) * mrSges(6,2) - m(5) * t235 - t115 * mrSges(5,1) - mrSges(5,3) * t262 + (t197 * t361 + t200 * t356) * t192 + t342 * (t191 * t261 - t193 * t291)) * g(3) + (m(4) * (-t107 * t199 - t108 * t196) * pkin(7) + t223 * t362 - t341) * qJD(3) + t304 * t101 + t355 * t118 + (t116 * t13 + t117 * t14 + (t276 * t92 + t307) * pkin(7) - t234 * t92 + t355 * t41 + (t101 + t81) * t40) * m(5) + (t199 * (-Ifges(4,2) * t196 + t317) + t211) * t270 / 0.2e1 + (t1 * t150 - t11 * t71 + t12 * t72 - t151 * t2) * mrSges(6,3) + (Ifges(6,1) * t151 + Ifges(6,4) * t150) * t333 + (-t107 * t276 - t108 * t277 + t350) * mrSges(4,3) + (-pkin(2) * t122 + t350 * pkin(7) - (t164 * t197 + (-t107 * t196 + t108 * t199) * t200) * t281) * m(4) - t344 * t277 + t354 * pkin(7) * t196 + t352 * t362 + t252 * t364 + t317 * t366 + t225 * t367 + (-t14 * mrSges(5,3) + t335) * t191 * t196 + (-t137 + t179) * mrSges(3,2) - t122 * t231 + (Ifges(6,5) * t151 + Ifges(6,6) * t150) * t327 + (-t233 - t267) * t172 + (t14 * mrSges(5,2) - t13 * mrSges(5,1) + pkin(7) * t138 + Ifges(4,4) * t366 + Ifges(4,2) * t367 + (t367 - t323) * Ifges(5,3) + (t373 - t326) * Ifges(5,6) + (t372 - t325) * Ifges(5,5)) * t199 + (-t13 * mrSges(5,3) + t331) * t193 * t196 + (-m(4) * t244 + mrSges(5,3) * t297 - t368 * (-qJ(4) * t297 - t144 * t322 + t244) + t342 * (-t144 * t292 - t145 * t193) + t361 * t145 + t336 * (-t144 * t288 + t145 * t191) + t337 * t144) * g(1) + (-m(4) * t245 + mrSges(5,3) * t298 - t368 * (-qJ(4) * t298 - t142 * t322 + t245) + t342 * (-t142 * t292 - t143 * t193) + t361 * t143 + t336 * (-t142 * t288 + t143 * t191) + t337 * t142) * g(2) + (Ifges(4,1) * t162 + Ifges(4,4) * t367 + t222 * t323 + t224 * t326 + t226 * t325) * t196 - t252 * t363 + (Ifges(6,4) * t151 + Ifges(6,2) * t150) * t332 + t353 * (pkin(7) * t276 - t234) + (t136 + t178) * mrSges(3,1) + t338 * t276 / 0.2e1 + Ifges(3,3) * qJDD(2) + t150 * t4 / 0.2e1 + t6 * (-mrSges(6,1) * t150 + mrSges(6,2) * t151) + t203 * t237 + qJDD(3) * (Ifges(4,5) * t196 + Ifges(4,6) * t199) - t69 * t252 / 0.2e1 + t154 * (Ifges(6,5) * t71 + Ifges(6,6) * t72 + Ifges(6,3) * t252) / 0.2e1 + t99 * (Ifges(6,4) * t71 + Ifges(6,2) * t72 + Ifges(6,6) * t252) / 0.2e1 + t159 * t258 + t229 * t307 + (Ifges(6,1) * t71 + Ifges(6,4) * t72 + Ifges(6,5) * t252) * t328 + t151 * t334; t315 * t325 + (mrSges(6,1) * t132 - mrSges(6,2) * t133 + t228 * t273) * t31 + t112 * t20 + t113 * t21 - pkin(3) * t52 - t36 * mrSges(4,2) + t37 * mrSges(4,1) + (t198 * t248 - t132 / 0.2e1) * t23 + (Ifges(6,5) * t133 + Ifges(6,6) * t132) * t324 + t357 * t56 + t358 * t55 + (-m(5) * t92 + t265 - t353) * t108 + ((qJ(4) * t6 + qJD(4) * t31) * m(6) - t40 * qJD(4) * m(5) + (t69 / 0.2e1 + t363 - t364 + Ifges(6,3) * t324 + Ifges(6,5) * t329 + Ifges(6,6) * t330) * t271 + t365 * qJ(4) + t6 * t227 + (Ifges(6,5) * t198 - Ifges(6,6) * t195) * t327 + t331 + (-Ifges(6,2) * t195 + t312) * t332 + (Ifges(6,1) * t198 - t313) * t333 + t369) * t191 + (Ifges(6,4) * t133 + Ifges(6,2) * t132) * t330 + t349 * mrSges(5,3) + t344 * t272 + (Ifges(5,2) * t326 + Ifges(5,6) * t323 - t335) * t193 + (Ifges(6,1) * t133 + Ifges(6,4) * t132) * t329 + (t147 * t370 + t148 * t339) * g(3) + (t339 * t86 + t370 * t85) * g(2) + (t339 * t88 + t370 * t87) * g(1) + (-t50 + t275) * t43 + t33 * t230 + (-t1 * t295 + t11 * t133 - t12 * t132 - t2 * t293 - t273 * t348) * mrSges(6,3) + (t195 * t248 - t133 / 0.2e1) * t24 + (-pkin(3) * t33 + t349 * qJ(4) + t375 * t41 - t40 * t57) * m(5) + t375 * t118 + (-t211 / 0.2e1 + t203 / 0.2e1) * t201 + (t154 * (-Ifges(6,5) * t195 - Ifges(6,6) * t198) - t213 * (-Ifges(6,1) * t195 - t312) + t99 * (-Ifges(6,2) * t198 - t313)) * t273 / 0.2e1 - (-Ifges(4,2) * t272 + t189 + t338) * t271 / 0.2e1 + (t1 * t113 + t357 * t11 + t112 * t2 + t358 * t12 - t31 * t50) * m(6) + Ifges(4,3) * qJDD(3) + (-t57 - t275) * t119 - t4 * t295 / 0.2e1 - Ifges(4,6) * t161 + Ifges(4,5) * t162 + t223 * t237 + (-t352 / 0.2e1 + t341) * qJD(2) + (t264 - t172) * t107 + t316 * t326 + t66 * t300 + t293 * t334; t195 * t21 + t198 * t20 + t304 * t158 + t219 * qJD(5) + (-t118 - t219) * t157 + t52 + (t1 * t195 + t154 * t348 - t158 * t31 + t198 * t2 + t210) * m(6) + (-t157 * t41 + t158 * t40 + t210 + t33) * m(5); -t31 * (-mrSges(6,1) * t213 + mrSges(6,2) * t99) + (Ifges(6,1) * t99 + t314) * t329 + t23 * t328 + (Ifges(6,5) * t99 + Ifges(6,6) * t213) * t324 - t11 * t55 + t12 * t56 - g(1) * ((-t195 * t45 + t198 * t87) * mrSges(6,1) + (-t195 * t87 - t198 * t45) * mrSges(6,2)) - g(2) * ((-t195 * t44 + t198 * t85) * mrSges(6,1) + (-t195 * t85 - t198 * t44) * mrSges(6,2)) - g(3) * (mrSges(6,1) * t46 - mrSges(6,2) * t47) + (t11 * t99 - t12 * t213) * mrSges(6,3) + t3 + (Ifges(6,2) * t213 + t24 + t95) * t330 + t340;];
tau = t5;
