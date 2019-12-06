% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:14
% EndTime: 2019-12-05 18:58:20
% DurationCPUTime: 3.86s
% Computational Cost: add. (9984->343), mult. (20155->431), div. (0->0), fcn. (18011->8), ass. (0->213)
t404 = Ifges(5,1) - Ifges(5,2);
t327 = Ifges(6,1) - Ifges(6,2);
t206 = sin(qJ(5));
t207 = sin(qJ(4));
t210 = cos(qJ(5));
t211 = cos(qJ(4));
t180 = -t206 * t207 + t210 * t211;
t394 = t180 / 0.2e1;
t212 = cos(qJ(3));
t204 = t207 ^ 2;
t205 = t211 ^ 2;
t284 = t204 + t205;
t364 = t284 * t212;
t208 = sin(qJ(3));
t335 = t208 * pkin(2);
t199 = pkin(8) + t335;
t325 = pkin(9) + t199;
t177 = t325 * t207;
t178 = t325 * t211;
t118 = -t177 * t206 + t178 * t210;
t263 = -t210 * t177 - t178 * t206;
t403 = -t118 * mrSges(6,1) - t263 * mrSges(6,2);
t347 = -pkin(9) - pkin(8);
t192 = t347 * t207;
t193 = t347 * t211;
t139 = t192 * t206 - t193 * t210;
t262 = t210 * t192 + t193 * t206;
t402 = -t139 * mrSges(6,1) - t262 * mrSges(6,2);
t209 = sin(qJ(2));
t194 = t208 * t209 * pkin(1);
t213 = cos(qJ(2));
t339 = pkin(1) * t213;
t201 = pkin(2) + t339;
t164 = t201 * t212 - t194;
t300 = t211 * mrSges(5,2);
t301 = t207 * mrSges(5,1);
t249 = t300 / 0.2e1 + t301 / 0.2e1;
t351 = m(6) * pkin(4);
t283 = t351 / 0.2e1;
t181 = -t206 * t211 - t210 * t207;
t86 = t181 * t164;
t388 = t86 / 0.2e1;
t396 = -mrSges(6,2) / 0.2e1;
t87 = t180 * t164;
t323 = mrSges(6,1) * t388 + t87 * t396;
t399 = -t249 * t164 + (t206 * t87 + t210 * t86) * t283 + t323;
t169 = t212 * t339 - t194;
t98 = t180 * t169;
t387 = -t98 / 0.2e1;
t97 = t181 * t169;
t322 = t97 * mrSges(6,1) / 0.2e1 + mrSges(6,2) * t387;
t398 = -t249 * t169 + (t206 * t98 + t210 * t97) * t283 + t322;
t334 = t212 * pkin(2);
t154 = t180 * t334;
t153 = t181 * t334;
t395 = t153 / 0.2e1;
t367 = mrSges(6,1) * t395 + t154 * t396;
t397 = -t249 * t334 + t367;
t285 = Ifges(6,5) * t180 + Ifges(6,6) * t181;
t31 = t285 + t403;
t392 = t31 * qJD(5);
t36 = t285 + t402;
t391 = t36 * qJD(5);
t287 = t209 * t212;
t165 = pkin(1) * t287 + t201 * t208;
t159 = pkin(8) + t165;
t326 = pkin(9) + t159;
t134 = t326 * t207;
t135 = t326 * t211;
t66 = -t134 * t206 + t135 * t210;
t385 = t66 * mrSges(6,1);
t390 = t285 - t385;
t172 = Ifges(6,4) * t180;
t315 = Ifges(6,4) * t181;
t358 = 0.2e1 * t172 * t394 + (-0.2e1 * t327 * t394 - t315) * t181;
t302 = t181 * mrSges(6,3);
t56 = t66 * t302;
t123 = -mrSges(6,1) * t181 + mrSges(6,2) * t180;
t158 = -pkin(3) - t164;
t336 = pkin(4) * t211;
t146 = t158 - t336;
t73 = t146 * t123;
t372 = t56 / 0.2e1 + t73 / 0.2e1;
t389 = t358 + t372;
t337 = pkin(4) * t207;
t386 = m(6) * t337;
t384 = t66 * mrSges(6,3);
t382 = t118 * mrSges(6,3);
t380 = t139 * mrSges(6,3);
t168 = (t208 * t213 + t287) * pkin(1);
t202 = -pkin(3) - t336;
t124 = -mrSges(6,1) * t180 - mrSges(6,2) * t181;
t186 = -mrSges(5,1) * t211 + mrSges(5,2) * t207;
t259 = t124 / 0.2e1 - mrSges(4,1) / 0.2e1 + t186 / 0.2e1;
t265 = t284 * t169;
t379 = -m(5) * (-pkin(3) * t168 + pkin(8) * t265) / 0.2e1 - m(6) * (t139 * t98 + t168 * t202 + t262 * t97) / 0.2e1 - t259 * t168;
t187 = t300 + t301;
t132 = t158 * t187;
t303 = t180 * mrSges(6,3);
t65 = -t210 * t134 - t206 * t135;
t55 = t65 * t303;
t376 = t132 / 0.2e1 + t55 / 0.2e1;
t200 = -pkin(3) - t334;
t166 = t200 * t187;
t80 = t263 * t303;
t375 = t166 / 0.2e1 + t80 / 0.2e1;
t185 = t202 - t334;
t105 = t185 * t123;
t81 = t118 * t302;
t371 = t81 / 0.2e1 + t105 / 0.2e1;
t122 = t124 * t337;
t316 = Ifges(5,4) * t211;
t355 = t316 * t211 + t122 + (-Ifges(5,4) * t207 + t211 * t404) * t207;
t352 = m(6) / 0.2e1;
t374 = 0.2e1 * t352;
t353 = m(5) / 0.2e1;
t373 = 0.2e1 * t353;
t254 = t65 * mrSges(6,2);
t230 = -t254 / 0.2e1;
t360 = t284 * mrSges(5,3);
t258 = -mrSges(4,2) + t360;
t281 = -mrSges(4,1) + t124 + t186;
t363 = t258 * t164 + (t87 * t180 + t86 * t181) * mrSges(6,3) + t281 * t165;
t106 = t262 * t303;
t338 = pkin(3) * t187;
t362 = -t106 / 0.2e1 + t338 / 0.2e1;
t107 = t139 * t302;
t114 = t202 * t123;
t361 = t107 / 0.2e1 + t114 / 0.2e1;
t359 = (t153 * t181 + t154 * t180) * mrSges(6,3);
t274 = t302 / 0.2e1;
t277 = -t303 / 0.2e1;
t357 = -t66 * t274 + t65 * t277 + t376;
t356 = -t118 * t274 + t263 * t277 + t375;
t354 = t258 * t169 + t281 * t168 + (t98 * t180 + t97 * t181) * mrSges(6,3) + (-mrSges(3,1) * t209 - mrSges(3,2) * t213) * pkin(1);
t349 = -mrSges(4,2) / 0.2e1;
t348 = -t66 / 0.2e1;
t346 = -t118 / 0.2e1;
t344 = -t139 / 0.2e1;
t343 = -t164 / 0.2e1;
t342 = t169 / 0.2e1;
t340 = m(6) * t208;
t324 = -t56 - t73;
t314 = pkin(4) * qJD(4);
t297 = -t105 - t81;
t225 = t204 * Ifges(5,4) + (-t404 * t207 - t316) * t211 - t122;
t257 = t327 * t181 - t172;
t13 = -t55 - t146 * t386 - t132 + (t315 + t384) * t181 + (t65 * mrSges(6,3) + t257) * t180 + t225 + t324;
t296 = t13 * qJD(1);
t266 = t284 * t164;
t16 = m(6) * (t146 * t165 + t65 * t86 + t66 * t87) + m(5) * (t158 * t165 + t159 * t266) + t363;
t293 = t16 * qJD(1);
t17 = m(6) * (t146 * t168 + t65 * t97 + t66 * t98) + m(5) * (t158 * t168 + t159 * t265) + m(4) * (-t164 * t168 + t165 * t169) + t354;
t292 = t17 * qJD(1);
t291 = t180 * t172;
t245 = t327 * t180 + t315;
t19 = -t291 + (t245 + t384) * t181 + t324;
t290 = t19 * qJD(1);
t286 = -t107 - t114;
t282 = t210 * t303;
t276 = t303 / 0.2e1;
t275 = -t302 / 0.2e1;
t271 = t395 + t388;
t270 = t154 / 0.2e1 + t87 / 0.2e1;
t269 = t343 + t342;
t267 = (t185 + t202) * t207;
t264 = t284 * t199;
t214 = -t66 * t275 + t65 * t276 - t355 - t376 - t389;
t219 = (t146 + t185) * t386;
t1 = -t219 / 0.2e1 + t263 * t276 - t118 * t275 + t214 - t371 - t375 + t398;
t18 = -t80 - t166 - t185 * t386 + (t315 + t382) * t181 + (mrSges(6,3) * t263 + t257) * t180 + t225 + t297;
t253 = -t1 * qJD(1) - t18 * qJD(2);
t21 = -t291 + (t245 + t382) * t181 + t297;
t234 = t371 + t389;
t223 = (t348 + t346) * t302 + t234;
t9 = t223 - t322;
t252 = t9 * qJD(1) - t21 * qJD(2);
t224 = t281 * t208 + t258 * t212;
t30 = m(6) * (t118 * t154 + t153 * t263) + t359 + (t185 * t340 + m(5) * (t364 * t199 + t200 * t208) + t224) * pkin(2);
t215 = (t259 * t208 + t212 * t349) * pkin(2) + t259 * t165 + (t165 * t200 + t164 * t264 + (t158 * t208 + t159 * t364) * pkin(2)) * t353 + (t118 * t87 + t146 * t335 + t153 * t65 + t154 * t66 + t165 * t185 + t263 * t86) * t352;
t4 = t269 * mrSges(4,2) + ((-t97 / 0.2e1 + t271) * t181 + (t387 + t270) * t180) * mrSges(6,3) + t215 + (t334 / 0.2e1 - t269) * t360 + t379;
t251 = t4 * qJD(1) + t30 * qJD(2);
t250 = t358 + t361;
t244 = m(6) * (t153 * t210 + t154 * t206);
t235 = t206 * pkin(4) * t302 + Ifges(5,5) * t211 - Ifges(5,6) * t207 + t285;
t20 = -t106 + t338 - t202 * t386 + (t315 + t380) * t181 + (mrSges(6,3) * t262 + t257) * t180 + t225 + t286;
t218 = (t146 + t202) * t386;
t5 = -t218 / 0.2e1 + t262 * t276 - t139 * t275 + t214 - t361 + t362 + t399;
t220 = -t139 * t274 + t262 * t277 + t250 + t355 - t362;
t216 = t220 + t356 + t371;
t7 = 0.2e1 * (-t244 / 0.4e1 + m(6) * t267 / 0.4e1) * pkin(4) + t216 - t397;
t233 = -t5 * qJD(1) + t7 * qJD(2) - t20 * qJD(3);
t222 = (t348 + t344) * t302 + t250 + t372;
t11 = t222 - t323;
t221 = (t346 + t344) * t302 + t250 + t371;
t14 = t221 - t367;
t25 = -t291 + (t245 + t380) * t181 + t286;
t232 = t11 * qJD(1) + t14 * qJD(2) - t25 * qJD(3);
t184 = (mrSges(6,1) * t206 + mrSges(6,2) * t210) * pkin(4);
t24 = t254 / 0.2e1 + (t348 + t66 / 0.2e1) * mrSges(6,1) + t230;
t33 = (t346 + t118 / 0.2e1) * mrSges(6,1);
t38 = (t344 + t139 / 0.2e1) * mrSges(6,1);
t226 = -qJD(1) * t24 - qJD(2) * t33 - qJD(3) * t38 + qJD(4) * t184;
t179 = t184 * qJD(5);
t22 = 0.2e1 * t230 + t390;
t15 = t221 + t367;
t12 = t222 + t323;
t10 = t223 + t322;
t8 = t216 + (t267 * t352 + t244 / 0.2e1) * pkin(4) + t397;
t6 = t218 / 0.2e1 + t220 + t357 + t372 + t399;
t3 = t97 * t274 + t98 * t276 + t169 * t349 + mrSges(4,2) * t343 + (t180 * t270 + t181 * t271) * mrSges(6,3) + t215 + t342 * t360 + (t164 + t334) * mrSges(5,3) * (t204 / 0.2e1 + t205 / 0.2e1) - t379;
t2 = t234 + t219 / 0.2e1 + t355 + t356 + t357 + t398;
t23 = [qJD(2) * t17 + qJD(3) * t16 - qJD(4) * t13 - qJD(5) * t19, t3 * qJD(3) + t2 * qJD(4) + t10 * qJD(5) + t292 + ((t118 * t98 + t168 * t185 + t263 * t97) * t374 + (t168 * t200 + t169 * t264) * t373 + m(4) * (-t168 * t212 + t169 * t208) * pkin(2) + t354) * qJD(2), t3 * qJD(2) + t6 * qJD(4) + t12 * qJD(5) + t293 + ((t139 * t87 + t165 * t202 + t262 * t86) * t374 + (-pkin(3) * t165 + pkin(8) * t266) * t373 + t363) * qJD(3), -t296 + t2 * qJD(2) + t6 * qJD(3) + (-pkin(4) * t282 - t254 - t385 + (t206 * t65 - t210 * t66) * t351 + t235 + t186 * t159) * qJD(4) + t22 * qJD(5), -t290 + t10 * qJD(2) + t12 * qJD(3) + t22 * qJD(4) + (-t254 + t390) * qJD(5); qJD(3) * t4 - qJD(4) * t1 + qJD(5) * t9 - t292, qJD(3) * t30 - qJD(4) * t18 - qJD(5) * t21, t8 * qJD(4) + t15 * qJD(5) + t251 + (m(6) * (t139 * t154 + t153 * t262) + t359 + (t202 * t340 + m(5) * (-pkin(3) * t208 + t364 * pkin(8)) + t224) * pkin(2)) * qJD(3), t8 * qJD(3) + (t186 * t199 + t235 + t403) * qJD(4) + t392 + (-t282 + m(6) * (-t118 * t210 + t206 * t263)) * t314 + t253, t15 * qJD(3) + t31 * qJD(4) + t252 + t392; -qJD(2) * t4 - qJD(4) * t5 + qJD(5) * t11 - t293, qJD(4) * t7 + qJD(5) * t14 - t251, -qJD(4) * t20 - qJD(5) * t25, (t186 * pkin(8) + t235 + t402) * qJD(4) + t391 + (-t282 + m(6) * (-t139 * t210 + t206 * t262)) * t314 + t233, t36 * qJD(4) + t232 + t391; qJD(2) * t1 + qJD(3) * t5 + qJD(5) * t24 + t296, -qJD(3) * t7 + qJD(5) * t33 - t253, qJD(5) * t38 - t233, -t179, -t179 - t226; -qJD(2) * t9 - qJD(3) * t11 - qJD(4) * t24 + t290, -qJD(3) * t14 - qJD(4) * t33 - t252, -qJD(4) * t38 - t232, t226, 0;];
Cq = t23;
