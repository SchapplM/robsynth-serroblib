% Calculate vector of inverse dynamics joint torques for
% S5RRRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:02
% EndTime: 2019-12-31 21:13:24
% DurationCPUTime: 11.42s
% Computational Cost: add. (8957->576), mult. (21020->783), div. (0->0), fcn. (15164->14), ass. (0->280)
t250 = qJ(2) + qJ(3);
t243 = pkin(9) + t250;
t230 = sin(t243);
t231 = cos(t243);
t253 = sin(qJ(5));
t364 = mrSges(6,2) * t253;
t434 = t230 * t364 + t231 * (m(6) * pkin(8) + mrSges(6,3));
t249 = qJD(2) + qJD(3);
t257 = cos(qJ(5));
t254 = sin(qJ(3));
t255 = sin(qJ(2));
t258 = cos(qJ(3));
t259 = cos(qJ(2));
t199 = -t254 * t255 + t258 * t259;
t184 = t199 * qJD(1);
t200 = t254 * t259 + t255 * t258;
t185 = t200 * qJD(1);
t251 = sin(pkin(9));
t252 = cos(pkin(9));
t281 = t184 * t251 + t185 * t252;
t123 = t249 * t257 - t253 * t281;
t124 = t249 * t253 + t257 * t281;
t411 = mrSges(5,1) * t249 + mrSges(6,1) * t123 - mrSges(6,2) * t124 - mrSges(5,3) * t281;
t261 = -pkin(7) - pkin(6);
t223 = t261 * t259;
t206 = qJD(1) * t223;
t186 = t254 * t206;
t222 = t261 * t255;
t205 = qJD(1) * t222;
t193 = qJD(2) * pkin(2) + t205;
t143 = t193 * t258 + t186;
t178 = t185 * qJ(4);
t117 = t143 - t178;
t108 = pkin(3) * t249 + t117;
t189 = t258 * t206;
t144 = t193 * t254 - t189;
t347 = qJ(4) * t184;
t118 = t144 + t347;
t344 = t118 * t251;
t61 = t108 * t252 - t344;
t58 = -pkin(4) * t249 - t61;
t433 = m(5) * t61 - m(6) * t58 + t411;
t111 = t252 * t118;
t62 = t108 * t251 + t111;
t59 = pkin(8) * t249 + t62;
t246 = t259 * pkin(2);
t236 = t246 + pkin(1);
t221 = t236 * qJD(1);
t156 = -pkin(3) * t184 + qJD(4) - t221;
t299 = t184 * t252 - t185 * t251;
t68 = -pkin(4) * t299 - pkin(8) * t281 + t156;
t20 = -t253 * t59 + t257 * t68;
t432 = t20 * mrSges(6,1);
t21 = t253 * t68 + t257 * t59;
t431 = t21 * mrSges(6,2);
t244 = sin(t250);
t245 = cos(t250);
t430 = mrSges(4,1) * t244 + mrSges(5,1) * t230 + mrSges(4,2) * t245 + mrSges(5,2) * t231;
t321 = qJD(1) * qJD(2);
t209 = qJDD(1) * t259 - t255 * t321;
t350 = t257 * mrSges(6,1);
t429 = t350 - t364;
t210 = qJDD(1) * t255 + t259 * t321;
t272 = t199 * qJD(3);
t126 = qJD(1) * t272 + t209 * t254 + t210 * t258;
t247 = qJDD(2) + qJDD(3);
t196 = t210 * pkin(6);
t155 = qJDD(2) * pkin(2) - pkin(7) * t210 - t196;
t195 = t209 * pkin(6);
t157 = pkin(7) * t209 + t195;
t83 = -qJD(3) * t144 + t155 * t258 - t157 * t254;
t39 = pkin(3) * t247 - qJ(4) * t126 - qJD(4) * t185 + t83;
t273 = t200 * qJD(3);
t127 = -qJD(1) * t273 + t209 * t258 - t210 * t254;
t324 = qJD(3) * t258;
t325 = qJD(3) * t254;
t82 = t155 * t254 + t157 * t258 + t193 * t324 + t206 * t325;
t44 = qJ(4) * t127 + qJD(4) * t184 + t82;
t15 = t251 * t39 + t252 * t44;
t12 = pkin(8) * t247 + t15;
t77 = -t126 * t251 + t127 * t252;
t78 = t126 * t252 + t127 * t251;
t346 = qJDD(1) * pkin(1);
t179 = -pkin(2) * t209 - t346;
t98 = -pkin(3) * t127 + qJDD(4) + t179;
t17 = -pkin(4) * t77 - pkin(8) * t78 + t98;
t2 = qJD(5) * t20 + t12 * t257 + t17 * t253;
t3 = -qJD(5) * t21 - t12 * t253 + t17 * t257;
t428 = t2 * t257 - t253 * t3;
t256 = sin(qJ(1));
t260 = cos(qJ(1));
t427 = g(1) * t260 + g(2) * t256;
t426 = -t245 * mrSges(4,1) - t231 * mrSges(5,1) + t244 * mrSges(4,2) + (mrSges(5,2) - mrSges(6,3)) * t230;
t45 = qJD(5) * t123 + t247 * t253 + t257 * t78;
t395 = t45 / 0.2e1;
t46 = -qJD(5) * t124 + t247 * t257 - t253 * t78;
t394 = t46 / 0.2e1;
t73 = qJDD(5) - t77;
t392 = t73 / 0.2e1;
t422 = m(5) + m(6);
t421 = t209 / 0.2e1;
t385 = t247 / 0.2e1;
t381 = t259 / 0.2e1;
t420 = -t281 / 0.2e1;
t419 = t281 / 0.2e1;
t418 = -t299 / 0.2e1;
t417 = Ifges(5,4) * t281;
t416 = Ifges(5,4) * t299;
t415 = Ifges(4,5) * t200;
t414 = Ifges(4,6) * t199;
t413 = t259 * Ifges(3,2);
t289 = mrSges(6,1) * t253 + mrSges(6,2) * t257;
t412 = t289 * t58;
t128 = -mrSges(5,2) * t249 + mrSges(5,3) * t299;
t132 = qJD(5) - t299;
t85 = -mrSges(6,2) * t132 + mrSges(6,3) * t123;
t86 = mrSges(6,1) * t132 - mrSges(6,3) * t124;
t282 = -t253 * t86 + t257 * t85;
t409 = -t128 - t282;
t159 = t222 * t254 - t223 * t258;
t295 = pkin(4) * t231 + pkin(8) * t230;
t152 = qJD(2) * t199 + t272;
t153 = -qJD(2) * t200 - t273;
t100 = t152 * t252 + t153 * t251;
t147 = t199 * t251 + t200 * t252;
t322 = qJD(5) * t257;
t276 = t100 * t253 + t147 * t322;
t327 = qJD(1) * t259;
t328 = qJD(1) * t255;
t371 = pkin(6) * t259;
t372 = pkin(6) * t255;
t407 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t328) * t371 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t327) * t372;
t406 = t195 * t259 + t196 * t255;
t404 = 0.2e1 * t385;
t403 = -t231 * t429 + t426;
t401 = -m(3) * pkin(6) + m(4) * t261 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t400 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t219 = -mrSges(3,1) * t259 + mrSges(3,2) * t255;
t399 = m(3) * pkin(1) + m(4) * t236 + mrSges(2,1) - t219 - t426;
t18 = mrSges(6,1) * t73 - mrSges(6,3) * t45;
t19 = -mrSges(6,2) * t73 + mrSges(6,3) * t46;
t323 = qJD(5) * t253;
t398 = m(6) * ((-t20 * t257 - t21 * t253) * qJD(5) + t428) - t86 * t322 - t85 * t323 + t257 * t19 - t253 * t18;
t396 = Ifges(6,1) * t395 + Ifges(6,4) * t394 + Ifges(6,5) * t392;
t391 = -t123 / 0.2e1;
t390 = -t124 / 0.2e1;
t389 = t124 / 0.2e1;
t388 = -t132 / 0.2e1;
t386 = t185 / 0.2e1;
t384 = -t249 / 0.2e1;
t382 = t257 / 0.2e1;
t380 = mrSges(5,3) * t62;
t379 = pkin(2) * t255;
t378 = pkin(2) * t258;
t377 = pkin(3) * t185;
t376 = pkin(3) * t244;
t234 = pkin(3) * t245;
t375 = pkin(3) * t251;
t374 = pkin(3) * t252;
t373 = pkin(4) * t230;
t366 = t61 * mrSges(5,3);
t363 = Ifges(3,4) * t255;
t362 = Ifges(3,4) * t259;
t361 = Ifges(6,4) * t253;
t360 = Ifges(6,4) * t257;
t359 = pkin(2) * qJD(3);
t358 = t124 * Ifges(6,4);
t357 = t143 * mrSges(4,3);
t356 = t185 * mrSges(4,3);
t355 = t185 * Ifges(4,4);
t343 = t299 * t253;
t342 = t299 * t257;
t341 = t147 * t253;
t340 = t147 * t257;
t335 = t251 * t254;
t334 = t252 * t254;
t333 = t253 * t256;
t332 = t253 * t260;
t331 = t256 * t257;
t330 = t257 * t260;
t151 = t205 * t258 + t186;
t235 = pkin(3) + t378;
t177 = pkin(2) * t334 + t235 * t251;
t329 = t234 + t246;
t326 = qJD(2) * t255;
t319 = Ifges(6,5) * t45 + Ifges(6,6) * t46 + Ifges(6,3) * t73;
t240 = pkin(2) * t326;
t315 = t230 * t350;
t239 = pkin(2) * t328;
t119 = Ifges(6,4) * t123;
t52 = Ifges(6,1) * t124 + t132 * Ifges(6,5) + t119;
t311 = t52 * t382;
t310 = t234 + t295;
t309 = qJD(2) * t261;
t306 = -t77 * mrSges(5,1) + mrSges(5,2) * t78;
t303 = -t323 / 0.2e1;
t140 = -pkin(3) * t153 + t240;
t302 = t321 / 0.2e1;
t150 = -t205 * t254 + t189;
t158 = t222 * t258 + t223 * t254;
t163 = -pkin(3) * t199 - t236;
t297 = t434 * t256;
t296 = t434 * t260;
t294 = -g(1) * t256 + g(2) * t260;
t293 = mrSges(3,1) * t255 + mrSges(3,2) * t259;
t288 = Ifges(6,1) * t257 - t361;
t287 = t363 + t413;
t286 = -Ifges(6,2) * t253 + t360;
t285 = Ifges(3,5) * t259 - Ifges(3,6) * t255;
t284 = Ifges(6,5) * t257 - Ifges(6,6) * t253;
t283 = -t20 * t253 + t21 * t257;
t14 = -t251 * t44 + t252 * t39;
t146 = -t199 * t252 + t200 * t251;
t87 = pkin(4) * t146 - pkin(8) * t147 + t163;
t135 = qJ(4) * t199 + t159;
t277 = -qJ(4) * t200 + t158;
t89 = t135 * t252 + t251 * t277;
t31 = -t253 * t89 + t257 * t87;
t32 = t253 * t87 + t257 * t89;
t176 = -pkin(2) * t335 + t235 * t252;
t279 = pkin(1) * t293;
t278 = t150 - t347;
t275 = -t100 * t257 + t147 * t323;
t274 = t255 * (Ifges(3,1) * t259 - t363);
t207 = t255 * t309;
t208 = t259 * t309;
t104 = t207 * t258 + t208 * t254 + t222 * t324 + t223 * t325;
t84 = pkin(4) * t281 - pkin(8) * t299 + t377;
t211 = -t376 - t379;
t266 = m(6) * (t211 - t373) - t315;
t265 = m(6) * (-t373 - t376) - t315;
t105 = -qJD(3) * t159 - t207 * t254 + t208 * t258;
t263 = -qJ(4) * t152 - qJD(4) * t200 + t105;
t11 = -pkin(4) * t247 - t14;
t133 = t184 * Ifges(4,2) + t249 * Ifges(4,6) + t355;
t180 = Ifges(4,4) * t184;
t134 = t185 * Ifges(4,1) + t249 * Ifges(4,5) + t180;
t50 = Ifges(6,5) * t124 + t123 * Ifges(6,6) + t132 * Ifges(6,3);
t51 = t123 * Ifges(6,2) + t132 * Ifges(6,6) + t358;
t8 = Ifges(6,4) * t45 + Ifges(6,2) * t46 + Ifges(6,6) * t73;
t91 = Ifges(5,2) * t299 + t249 * Ifges(5,6) + t417;
t92 = Ifges(5,1) * t281 + t249 * Ifges(5,5) + t416;
t262 = t221 * (mrSges(4,1) * t185 + mrSges(4,2) * t184) + (-t417 + t50) * t420 + (-t156 * mrSges(5,2) + Ifges(5,1) * t420 + Ifges(5,5) * t384 + t284 * t388 + t286 * t391 + t288 * t390 + t366 - t412) * t299 + t144 * t356 + t184 * t357 + (Ifges(4,3) + Ifges(5,3)) * t247 + (t311 + t412) * qJD(5) - t15 * mrSges(5,2) + t14 * mrSges(5,1) + (t416 + t92) * t418 + (t303 + t343 / 0.2e1) * t51 + Ifges(5,6) * t77 + Ifges(5,5) * t78 - t82 * mrSges(4,2) + t83 * mrSges(4,1) + t8 * t382 + (Ifges(4,5) * t184 - Ifges(4,6) * t185) * t384 + t133 * t386 + (Ifges(6,5) * t253 + Ifges(6,6) * t257) * t392 + (Ifges(6,2) * t257 + t361) * t394 + (Ifges(6,1) * t253 + t360) * t395 + t253 * t396 + t91 * t419 + Ifges(4,5) * t126 + Ifges(4,6) * t127 + ((-t323 + t343) * t21 + (-t322 + t342) * t20 + t428) * mrSges(6,3) - t11 * t429 + (-t156 * mrSges(5,1) + Ifges(6,5) * t390 - Ifges(5,2) * t418 - Ifges(5,6) * t384 + Ifges(6,6) * t391 + Ifges(6,3) * t388 + t380 + t431 - t432) * t281 - t52 * t342 / 0.2e1 - t185 * (Ifges(4,1) * t184 - t355) / 0.2e1 + (t123 * t286 + t124 * t288 + t132 * t284) * qJD(5) / 0.2e1 - (-Ifges(4,2) * t185 + t134 + t180) * t184 / 0.2e1;
t248 = -qJ(4) + t261;
t238 = Ifges(3,4) * t327;
t233 = -pkin(4) - t374;
t204 = pkin(1) + t329;
t183 = Ifges(3,1) * t328 + Ifges(3,5) * qJD(2) + t238;
t182 = Ifges(3,6) * qJD(2) + qJD(1) * t287;
t171 = -pkin(4) - t176;
t170 = t231 * t330 + t333;
t169 = -t231 * t332 + t331;
t168 = -t231 * t331 + t332;
t167 = t231 * t333 + t330;
t162 = mrSges(4,1) * t249 - t356;
t161 = -mrSges(4,2) * t249 + mrSges(4,3) * t184;
t160 = t239 + t377;
t142 = -mrSges(4,1) * t184 + mrSges(4,2) * t185;
t125 = -t178 + t151;
t110 = -mrSges(4,2) * t247 + mrSges(4,3) * t127;
t109 = mrSges(4,1) * t247 - mrSges(4,3) * t126;
t99 = t152 * t251 - t153 * t252;
t97 = -mrSges(5,1) * t299 + mrSges(5,2) * t281;
t79 = t239 + t84;
t76 = t125 * t252 + t251 * t278;
t74 = qJ(4) * t153 + qJD(4) * t199 + t104;
t66 = t117 * t252 - t344;
t65 = t117 * t251 + t111;
t64 = mrSges(5,1) * t247 - mrSges(5,3) * t78;
t63 = -mrSges(5,2) * t247 + mrSges(5,3) * t77;
t37 = pkin(4) * t99 - pkin(8) * t100 + t140;
t27 = t253 * t79 + t257 * t76;
t26 = -t253 * t76 + t257 * t79;
t25 = t253 * t84 + t257 * t66;
t24 = -t253 * t66 + t257 * t84;
t23 = t251 * t263 + t252 * t74;
t16 = -mrSges(6,1) * t46 + mrSges(6,2) * t45;
t5 = -qJD(5) * t32 - t23 * t253 + t257 * t37;
t4 = qJD(5) * t31 + t23 * t257 + t253 * t37;
t1 = [-t433 * (t251 * t74 - t252 * t263) + (-m(5) * t14 + m(6) * t11 + t16 - t64) * (t135 * t251 - t252 * t277) + (Ifges(4,5) * t152 + Ifges(5,5) * t100 + Ifges(4,6) * t153 - Ifges(5,6) * t99) * t249 / 0.2e1 + t299 * (Ifges(5,4) * t100 - Ifges(5,2) * t99) / 0.2e1 + (-mrSges(4,2) * t236 + Ifges(4,1) * t200 + Ifges(4,4) * t199) * t126 + (t259 * t362 + t274) * t302 + (Ifges(3,4) * t210 + Ifges(3,2) * t209) * t381 + t99 * t432 + (mrSges(4,1) * t236 + Ifges(4,4) * t200 + Ifges(4,2) * t199) * t127 + t100 * t311 + (t144 * t153 + t199 * t82 - t200 * t83) * mrSges(4,3) + t210 * t362 / 0.2e1 + (t98 * mrSges(5,2) - t14 * mrSges(5,3) + Ifges(5,1) * t78 + Ifges(5,4) * t77 + Ifges(5,5) * t404 + t11 * t289 + t284 * t392 + t286 * t394 + t288 * t395 + t52 * t303) * t147 + (-t15 * mrSges(5,3) + Ifges(6,3) * t392 + Ifges(6,6) * t394 + Ifges(6,5) * t395 - Ifges(5,2) * t77 - Ifges(5,4) * t78 + t98 * mrSges(5,1) + t319 / 0.2e1 - t404 * Ifges(5,6) + t400) * t146 + (t209 * t371 + t210 * t372 + t406) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t406) + (t183 * t381 + t285 * qJD(2) / 0.2e1 - t407) * qJD(2) - t276 * t51 / 0.2e1 + (-t2 * t341 + t20 * t275 - t21 * t276 - t3 * t340) * mrSges(6,3) + m(5) * (t140 * t156 + t15 * t89 + t163 * t98 + t23 * t62) + m(6) * (t2 * t32 + t20 * t5 + t21 * t4 + t3 * t31) + t142 * t240 + (t414 / 0.2e1 + t415 / 0.2e1) * t247 + (t414 + t415) * t385 + t31 * t18 + t32 * t19 + (-mrSges(3,1) * t372 - mrSges(3,2) * t371 + 0.2e1 * Ifges(3,6) * t381) * qJDD(2) + (Ifges(3,1) * t210 + Ifges(3,4) * t421 + Ifges(3,5) * qJDD(2) - t302 * t413) * t255 + t58 * (mrSges(6,1) * t276 - mrSges(6,2) * t275) + t132 * (-Ifges(6,5) * t275 - Ifges(6,6) * t276 + Ifges(6,3) * t99) / 0.2e1 + (-t168 * mrSges(6,1) - t167 * mrSges(6,2) + (t248 * t422 + t401) * t260 + (m(5) * t204 - m(6) * (-t204 - t295) + t399) * t256) * g(1) + t123 * (-Ifges(6,4) * t275 - Ifges(6,2) * t276 + Ifges(6,6) * t99) / 0.2e1 + (-t170 * mrSges(6,1) - t169 * mrSges(6,2) - t422 * (t204 * t260 - t248 * t256) + t401 * t256 + (-m(6) * t295 - t399) * t260) * g(2) + t4 * t85 + t5 * t86 + t89 * t63 + t99 * t50 / 0.2e1 - t99 * t91 / 0.2e1 + t100 * t92 / 0.2e1 + (Ifges(4,1) * t152 + Ifges(4,4) * t153) * t386 + (-Ifges(6,1) * t275 - Ifges(6,4) * t276 + Ifges(6,5) * t99) * t389 + t340 * t396 + (Ifges(5,1) * t100 - Ifges(5,4) * t99) * t419 + t287 * t421 + t23 * t128 + t140 * t97 + t152 * t134 / 0.2e1 + t153 * t133 / 0.2e1 + t156 * (mrSges(5,1) * t99 + mrSges(5,2) * t100) + t158 * t109 + t159 * t110 + t104 * t161 + t105 * t162 + t184 * (Ifges(4,4) * t152 + Ifges(4,2) * t153) / 0.2e1 + t179 * (-mrSges(4,1) * t199 + mrSges(4,2) * t200) - pkin(1) * (-mrSges(3,1) * t209 + mrSges(3,2) * t210) - t221 * (-mrSges(4,1) * t153 + mrSges(4,2) * t152) + t163 * t306 - t99 * t431 - t279 * t321 - t182 * t326 / 0.2e1 - t8 * t341 / 0.2e1 - t219 * t346 - t152 * t357 + Ifges(2,3) * qJDD(1) - t100 * t366 + m(4) * (t104 * t144 + t105 * t143 + t158 * t83 + t159 * t82 - t179 * t236 - t221 * t240) - t99 * t380; t433 * (t125 * t251 - t252 * t278 - (t251 * t258 + t334) * t359) - (-Ifges(3,2) * t328 + t183 + t238) * t327 / 0.2e1 + (m(5) * t62 + m(6) * t283 - t409) * (t252 * t258 - t335) * t359 + (-t162 * t325 + t161 * t324 + t254 * t110 + m(4) * (t254 * t82 + t258 * t83 + (-t143 * t254 + t144 * t258) * qJD(3))) * pkin(2) + (t219 - m(6) * (t246 + t310) - m(5) * t329 - m(4) * t246 + t403) * g(3) + (t14 * t176 + t15 * t177 - t156 * t160 - t62 * t76) * m(5) + (t11 * t171 - t20 * t26 - t21 * t27) * m(6) + t398 * (pkin(8) + t177) + t262 + (t407 + (-t274 / 0.2e1 + t279) * qJD(1)) * qJD(1) - t27 * t85 - t26 * t86 - g(1) * (t260 * t266 + t296) - g(2) * (t256 * t266 + t297) + t109 * t378 - t76 * t128 - t160 * t97 - t151 * t161 - t150 * t162 + t171 * t16 + t176 * t64 + t177 * t63 - t195 * mrSges(3,2) - t196 * mrSges(3,1) + Ifges(3,3) * qJDD(2) + Ifges(3,6) * t209 + Ifges(3,5) * t210 + t427 * (m(4) * t379 - m(5) * t211 + t293 + t430) - t142 * t239 - t285 * t321 / 0.2e1 + t182 * t328 / 0.2e1 - m(4) * (t143 * t150 + t144 * t151 - t221 * t239); t64 * t374 + t63 * t375 + t262 - t25 * t85 - t24 * t86 - g(1) * (t260 * t265 + t296) - g(2) * (t256 * t265 + t297) - t66 * t128 - t143 * t161 + t144 * t162 + t233 * t16 - t97 * t377 + t411 * t65 + (t11 * t233 - t20 * t24 - t21 * t25 - t58 * t65) * m(6) + ((t14 * t252 + t15 * t251) * pkin(3) - t156 * t377 + t61 * t65 - t62 * t66) * m(5) + (-m(5) * t234 - m(6) * t310 + t403) * g(3) + t398 * (pkin(8) + t375) + (m(5) * t376 + t430) * t427; t257 * t18 + t253 * t19 + t411 * t281 + t282 * qJD(5) + t409 * t299 + t306 + (t132 * t283 + t2 * t253 + t257 * t3 - t281 * t58 + t294) * m(6) + (t281 * t61 - t299 * t62 + t294 + t98) * m(5); -t58 * (mrSges(6,1) * t124 + mrSges(6,2) * t123) + (Ifges(6,1) * t123 - t358) * t390 + t51 * t389 + (Ifges(6,5) * t123 - Ifges(6,6) * t124) * t388 - t20 * t85 + t21 * t86 - g(1) * (mrSges(6,1) * t169 - mrSges(6,2) * t170) - g(2) * (-mrSges(6,1) * t167 + mrSges(6,2) * t168) + g(3) * t289 * t230 + (t123 * t20 + t124 * t21) * mrSges(6,3) + t319 + (-Ifges(6,2) * t124 + t119 + t52) * t391 + t400;];
tau = t1;
