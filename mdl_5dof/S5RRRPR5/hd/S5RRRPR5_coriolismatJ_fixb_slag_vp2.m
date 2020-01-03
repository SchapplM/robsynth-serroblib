% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:06
% EndTime: 2019-12-31 21:13:16
% DurationCPUTime: 5.67s
% Computational Cost: add. (14485->343), mult. (28449->475), div. (0->0), fcn. (32617->8), ass. (0->220)
t230 = sin(qJ(3));
t231 = sin(qJ(2));
t233 = cos(qJ(3));
t234 = cos(qJ(2));
t205 = -t230 * t231 + t233 * t234;
t221 = -pkin(2) * t234 - pkin(1);
t182 = -t205 * pkin(3) + t221;
t438 = m(5) * t182;
t385 = -pkin(7) - pkin(6);
t301 = t385 * t231;
t399 = t385 * t234;
t407 = t230 * t399 + t233 * t301;
t414 = t407 * mrSges(4,2);
t175 = -t230 * t301 + t233 * t399;
t419 = t175 * mrSges(4,1);
t228 = sin(pkin(9));
t245 = t205 * qJ(4) - t175;
t329 = cos(pkin(9));
t206 = -t230 * t234 - t233 * t231;
t413 = t206 * qJ(4) + t407;
t423 = -t228 * t245 + t329 * t413;
t432 = t423 * mrSges(5,2);
t437 = -t414 / 0.2e1 + t419 / 0.2e1 - t432 / 0.2e1;
t229 = sin(qJ(5));
t232 = cos(qJ(5));
t270 = mrSges(6,1) * t232 - mrSges(6,2) * t229;
t422 = t228 * t413 + t329 * t245;
t426 = t422 * t270;
t433 = t422 * mrSges(5,1);
t436 = -t426 / 0.2e1 - t433 / 0.2e1;
t167 = t329 * t205 + t206 * t228;
t253 = t228 * t205 - t329 * t206;
t267 = Ifges(6,5) * t229 + Ifges(6,6) * t232;
t373 = t232 / 0.2e1;
t375 = t229 / 0.2e1;
t383 = t253 / 0.2e1;
t357 = Ifges(6,2) * t232;
t359 = Ifges(6,4) * t229;
t211 = t357 + t359;
t376 = -t229 / 0.2e1;
t201 = t211 * t376;
t223 = Ifges(6,4) * t232;
t361 = Ifges(6,1) * t229;
t213 = t223 + t361;
t408 = t213 * t373 + t201;
t395 = -Ifges(6,2) * t229 + t223;
t68 = Ifges(6,6) * t253 + t167 * t395;
t214 = Ifges(6,1) * t232 - t359;
t70 = Ifges(6,5) * t253 + t214 * t167;
t248 = Ifges(4,5) * t205 + Ifges(4,6) * t206 - Ifges(5,6) * t253 + t267 * t383 + t68 * t373 + t70 * t375 + (Ifges(5,5) + t408) * t167;
t435 = t248 - t414 + t419 - t426 - t432 - t433;
t340 = t232 * mrSges(6,2);
t342 = t229 * mrSges(6,1);
t210 = t340 + t342;
t117 = t210 * t167;
t118 = t210 * t253;
t325 = t167 * t229;
t119 = -mrSges(6,2) * t253 - mrSges(6,3) * t325;
t324 = t167 * t232;
t121 = mrSges(6,1) * t253 - mrSges(6,3) * t324;
t162 = t167 * mrSges(5,2);
t412 = Ifges(6,5) * t167;
t71 = t214 * t253 - t412;
t338 = t232 * t71;
t411 = Ifges(6,6) * t167;
t69 = t253 * t395 - t411;
t341 = t229 * t69;
t353 = t253 * mrSges(5,1);
t93 = -pkin(4) * t167 - pkin(8) * t253 + t182;
t39 = -t229 * t422 + t232 * t93;
t222 = Ifges(6,5) * t232;
t396 = -Ifges(6,6) * t229 + t222;
t40 = t229 * t93 + t232 * t422;
t405 = -t167 / 0.2e1;
t434 = t422 * t118 + t182 * (t162 + t353) + t221 * (-mrSges(4,1) * t206 + mrSges(4,2) * t205) + t39 * t121 + t40 * t119 - t117 * t423 + (t338 / 0.2e1 - t341 / 0.2e1 + t396 * t405 + Ifges(5,4) * t167) * t167 + (-Ifges(5,4) * t253 + t70 * t373 + t68 * t376 + t383 * t396 + (-Ifges(5,2) - Ifges(6,3) + Ifges(5,1)) * t167) * t253;
t220 = pkin(2) * t233 + pkin(3);
t316 = t228 * t230;
t189 = -pkin(2) * t316 + t329 * t220;
t186 = -pkin(4) - t189;
t431 = t186 * t422;
t280 = t329 * t230;
t190 = pkin(2) * t280 + t228 * t220;
t430 = t190 * t423;
t292 = t329 * pkin(3);
t217 = -t292 - pkin(4);
t429 = t217 * t422;
t428 = t229 * t423;
t427 = t232 * t423;
t367 = t422 * t423;
t425 = t228 * t423 - t329 * t422;
t410 = t167 * mrSges(5,3);
t226 = t229 ^ 2;
t227 = t232 ^ 2;
t309 = t226 + t227;
t409 = t309 * mrSges(6,3) - mrSges(5,2);
t398 = t119 * t375 + t121 * t373;
t404 = mrSges(5,3) * t253;
t402 = -mrSges(5,1) - t270;
t355 = Ifges(6,3) * t253;
t382 = -t167 / 0.4e1;
t394 = t222 * t382 - t355 / 0.2e1;
t371 = t206 * pkin(3);
t101 = pkin(4) * t253 - pkin(8) * t167 - t371;
t45 = t101 * t232 - t428;
t46 = t101 * t229 + t427;
t264 = -t45 * t229 + t46 * t232;
t393 = (-mrSges(4,1) * t230 - mrSges(4,2) * t233) * pkin(2);
t392 = -Ifges(6,6) * t325 / 0.2e1 + Ifges(6,5) * t324 / 0.2e1 + t355 / 0.2e1 + t396 * t382 + t338 / 0.4e1 - t341 / 0.4e1;
t362 = mrSges(6,3) * t253;
t305 = t229 * t362;
t120 = mrSges(6,2) * t167 - t305;
t122 = -t167 * mrSges(6,1) - t232 * t362;
t187 = pkin(8) + t190;
t195 = (t329 * t233 - t316) * pkin(2);
t339 = t232 * t40;
t266 = -t229 * t39 + t339;
t319 = t190 * t253;
t320 = t189 * t167;
t322 = t186 * t117;
t194 = (t228 * t233 + t280) * pkin(2);
t344 = t194 * t423;
t379 = -t195 / 0.2e1;
t380 = t194 / 0.2e1;
t381 = -t187 / 0.2e1;
t386 = t46 / 0.2e1;
t387 = -t45 / 0.2e1;
t389 = m(6) / 0.2e1;
t390 = m(5) / 0.2e1;
t391 = (t430 - t344 + (-t189 + t195) * t422) * t390 + (t264 * t187 + t266 * t195 - t344 + t431) * t389 + t322 / 0.2e1 + t118 * t380 + (mrSges(6,3) * t387 + t121 * t381 + t122 * t379) * t229 + (mrSges(6,3) * t386 + t195 * t120 / 0.2e1 + t187 * t119 / 0.2e1) * t232 + (-t167 * t379 + t253 * t380 - t320 / 0.2e1 - t319 / 0.2e1) * mrSges(5,3) + t437;
t388 = m(5) * pkin(3);
t116 = t270 * t253;
t384 = t116 / 0.2e1;
t378 = -t270 / 0.2e1;
t372 = pkin(3) * t228;
t216 = pkin(8) + t372;
t377 = -t216 / 0.2e1;
t374 = -t232 / 0.2e1;
t225 = t231 * pkin(2);
t360 = Ifges(4,4) * t206;
t125 = -mrSges(5,1) * t167 + mrSges(5,2) * t253;
t185 = t225 - t371;
t271 = Ifges(4,4) * t205 + (-Ifges(4,1) + Ifges(4,2)) * t206;
t94 = t101 + t225;
t41 = t232 * t94 - t428;
t42 = t229 * t94 + t427;
t1 = m(4) * t221 * t225 + m(6) * (t39 * t41 + t40 * t42 - t367) + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t231) * t231 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t234 + (Ifges(3,1) - Ifges(3,2)) * t231) * t234 + (-mrSges(4,1) * t225 + t271) * t205 + (-mrSges(4,2) * t225 - t360) * t206 + t42 * t120 + t41 * t122 + (t125 + t438) * t185 + t434;
t354 = t1 * qJD(1);
t352 = t253 * t423;
t349 = t167 * mrSges(6,3);
t3 = -t371 * t438 + m(6) * (t39 * t45 + t40 * t46 - t367) + t271 * t205 + (-pkin(3) * t125 - t360) * t206 + t46 * t120 + t45 * t122 + t434;
t337 = t3 * qJD(1);
t336 = t41 * t229;
t335 = t42 * t232;
t5 = t40 * t122 + t423 * t116 + (t71 * t375 + t69 * t373 + mrSges(6,3) * t339 + t267 * t405 + (-t213 * t374 + t201) * t253) * t253 + (-t120 - t305) * t39;
t332 = t5 * qJD(1);
t311 = t232 * t120;
t313 = t229 * t122;
t10 = (t118 + t404) * t253 + (t311 - t313 + t410) * t167 + m(6) * (t266 * t167 - t352) + m(5) * (t167 * t422 - t352);
t328 = qJD(1) * t10;
t278 = t309 * t187;
t281 = -t227 / 0.2e1 - t226 / 0.2e1;
t241 = -t281 * t349 + (t167 * t190 - t189 * t253) * t390 + (t167 * t278 + t186 * t253) * t389;
t247 = t185 * t390 + (t229 * t42 + t232 * t41) * t389;
t252 = -t162 - t398;
t12 = (t378 - mrSges(5,1)) * t253 + t241 - t247 + t252;
t327 = t12 * qJD(1);
t277 = t309 * t216;
t291 = t253 * t378;
t240 = (t167 * t277 + t217 * t253) * t389 + t291 + (t167 * t228 - t253 * t329) * t388 / 0.2e1 + t309 * t349 / 0.2e1;
t308 = -t388 / 0.2e1;
t244 = (t229 * t46 + t232 * t45) * t389 + t206 * t308;
t14 = t240 - t244 + t252 - t353;
t326 = t14 * qJD(1);
t258 = -t340 / 0.2e1 - t342 / 0.2e1;
t251 = t258 * t167;
t255 = t313 / 0.2e1 - t311 / 0.2e1;
t17 = t251 + t255;
t323 = t17 * qJD(1);
t321 = t186 * t210;
t318 = t217 * t117;
t317 = t217 * t210;
t314 = t229 * t121;
t312 = t232 * t119;
t307 = mrSges(6,3) * t336;
t306 = mrSges(6,3) * t335;
t299 = t216 * t314;
t298 = t216 * t312;
t293 = -t423 * t210 / 0.2e1;
t276 = t372 * t404;
t274 = t214 * t375 + t373 * t395 + t408;
t273 = t71 / 0.4e1 - t412 / 0.2e1;
t272 = mrSges(6,3) * t281;
t265 = t335 - t336;
t263 = t292 * t410;
t26 = t402 * t194 + t393 + t409 * t195 + m(5) * (-t189 * t194 + t190 * t195) + m(6) * (t186 * t194 + t195 * t278);
t235 = -m(6) * (t265 * t216 + t429) / 0.2e1 - t318 / 0.2e1 + t425 * t308 + t299 / 0.2e1 - t298 / 0.2e1 + t307 / 0.2e1 - t306 / 0.2e1 + t276 / 0.2e1 + t263 / 0.2e1 - t436 - t437;
t4 = t235 + t391 + t436;
t262 = t4 * qJD(1) + t26 * qJD(2);
t108 = t274 + t321;
t242 = (t214 / 0.4e1 - t211 / 0.4e1 - t357 / 0.4e1) * t232 + (-t213 / 0.4e1 - t395 / 0.4e1 - t223 / 0.2e1 - t361 / 0.4e1) * t229;
t239 = (t187 * t272 + t242) * t253 + t186 * t384 + t293;
t254 = 0.3e1 / 0.4e1 * t411 - t69 / 0.4e1;
t260 = -t41 * mrSges(6,1) / 0.2e1 + t42 * mrSges(6,2) / 0.2e1;
t7 = (t120 * t381 + t254) * t229 + (t122 * t381 + t273) * t232 + t239 + t260 + t394;
t261 = t7 * qJD(1) + t108 * qJD(2);
t259 = mrSges(6,1) * t387 + mrSges(6,2) * t386;
t256 = t120 * t376 + t122 * t374;
t123 = t274 + t317;
t57 = (-t186 / 0.2e1 - t217 / 0.2e1) * t210 + (mrSges(6,2) * t379 - t213 / 0.2e1 - t395 / 0.2e1) * t232 + (mrSges(6,1) * t379 - t214 / 0.2e1 + t211 / 0.2e1) * t229;
t238 = (t216 * t272 + t242) * t253 + t217 * t384 + t293;
t9 = (t120 * t377 + t254) * t229 + (t122 * t377 + t273) * t232 + t238 + t259 + t394;
t250 = t9 * qJD(1) - t57 * qJD(2) + t123 * qJD(3);
t58 = t321 / 0.2e1 + t317 / 0.2e1 + t258 * t195 + t274;
t18 = t251 - t255;
t15 = t240 + t244 + t398;
t13 = t291 + t241 + t247 + t398;
t8 = t256 * t216 + t238 - t259 + t392;
t6 = t256 * t187 + t239 - t260 + t392;
t2 = -t235 + t248 + (-mrSges(5,1) / 0.2e1 + t378) * t422 + t391;
t11 = [qJD(2) * t1 + qJD(3) * t3 + qJD(4) * t10 - qJD(5) * t5, t2 * qJD(3) + t13 * qJD(4) + t6 * qJD(5) + t354 + (Ifges(3,5) * t234 - Ifges(3,6) * t231 + t187 * t312 - t187 * t314 + t306 - t307 + t322 + 0.2e1 * (t265 * t187 + t431) * t389 + 0.2e1 * (-t189 * t422 + t430) * t390 + (m(4) * (t175 * t233 + t230 * t407) + (-t205 * t233 + t206 * t230) * mrSges(4,3)) * pkin(2) + (-mrSges(3,1) * t234 + mrSges(3,2) * t231) * pkin(6) + (-t319 - t320) * mrSges(5,3) + t435) * qJD(2), t337 + t2 * qJD(2) + (t298 - t299 + m(6) * (t216 * t264 + t429) + t425 * t388 - t263 - t276 + t318 + t264 * mrSges(6,3) + t435) * qJD(3) + t15 * qJD(4) + t8 * qJD(5), qJD(2) * t13 + qJD(3) * t15 + qJD(5) * t18 + t328, -t332 + t6 * qJD(2) + t8 * qJD(3) + t18 * qJD(4) + (-mrSges(6,1) * t40 - mrSges(6,2) * t39 - t253 * t267) * qJD(5); qJD(3) * t4 + qJD(4) * t12 + qJD(5) * t7 - t354, qJD(3) * t26 + qJD(5) * t108, ((m(6) * t217 - t329 * t388 + t402) * t194 + t393 + (m(6) * t277 + t228 * t388 + t409) * t195) * qJD(3) + t58 * qJD(5) + t262, t327, t58 * qJD(3) + (-t187 * t270 + t396) * qJD(5) + t261; -qJD(2) * t4 + qJD(4) * t14 + qJD(5) * t9 - t337, -qJD(5) * t57 - t262, t123 * qJD(5), t326, (-t216 * t270 + t396) * qJD(5) + t250; -qJD(2) * t12 - qJD(3) * t14 - qJD(5) * t17 - t328, -t327, -t326, 0, -qJD(5) * t210 - t323; -qJD(2) * t7 - qJD(3) * t9 + qJD(4) * t17 + t332, qJD(3) * t57 - t261, -t250, t323, 0;];
Cq = t11;
