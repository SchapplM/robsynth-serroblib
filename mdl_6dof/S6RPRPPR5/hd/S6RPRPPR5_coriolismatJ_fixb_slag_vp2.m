% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:09
% EndTime: 2019-03-09 02:50:17
% DurationCPUTime: 4.97s
% Computational Cost: add. (13985->410), mult. (27030->546), div. (0->0), fcn. (30870->8), ass. (0->227)
t232 = sin(pkin(9));
t234 = cos(pkin(9));
t236 = sin(qJ(3));
t355 = cos(qJ(3));
t216 = t232 * t355 + t236 * t234;
t231 = sin(pkin(10));
t233 = cos(pkin(10));
t237 = cos(qJ(6));
t354 = sin(qJ(6));
t262 = t237 * t231 + t233 * t354;
t254 = t262 * t216;
t263 = t231 * t354 - t237 * t233;
t383 = t263 * t216;
t396 = m(7) * (-t263 * t254 + t262 * t383);
t308 = t263 * t383;
t312 = t262 * t254;
t395 = t308 + t312;
t372 = m(7) / 0.2e1;
t363 = -t262 / 0.2e1;
t361 = t263 / 0.2e1;
t213 = t232 * t236 - t234 * t355;
t151 = t262 * t213;
t307 = t263 * t151;
t148 = t263 * t213;
t314 = t262 * t148;
t294 = (t307 - t314) * t372;
t350 = pkin(3) + qJ(5);
t393 = t216 * t350;
t348 = -pkin(8) - t350;
t218 = t348 * t231;
t219 = t348 * t233;
t172 = -t218 * t354 + t237 * t219;
t173 = t237 * t218 + t219 * t354;
t284 = -pkin(2) * t234 - pkin(1);
t268 = -qJ(4) * t216 + t284;
t120 = t213 * t350 + t268;
t349 = pkin(7) + qJ(2);
t220 = t349 * t232;
t221 = t349 * t234;
t184 = t355 * t220 + t221 * t236;
t141 = pkin(4) * t216 + t184;
t128 = t233 * t141;
t59 = -t120 * t231 + t128;
t60 = t233 * t120 + t231 * t141;
t275 = t231 * t60 + t233 * t59;
t48 = pkin(5) * t216 + t128 + (-pkin(8) * t213 - t120) * t231;
t311 = t213 * t233;
t50 = pkin(8) * t311 + t60;
t29 = t237 * t48 - t354 * t50;
t30 = t237 * t50 + t354 * t48;
t328 = t231 * mrSges(6,3);
t168 = t216 * mrSges(6,1) - t213 * t328;
t304 = t233 * t168;
t170 = -t216 * mrSges(6,2) + mrSges(6,3) * t311;
t305 = t231 * t170;
t103 = -mrSges(7,2) * t216 - t148 * mrSges(7,3);
t105 = mrSges(7,1) * t216 - t151 * mrSges(7,3);
t346 = t103 * t363 + t105 * t361;
t373 = -m(7) / 0.2e1;
t374 = m(6) / 0.2e1;
t392 = t275 * t374 + (-t173 * t148 - t172 * t151 - t262 * t30 + t263 * t29) * t373 + t305 / 0.2e1 + t304 / 0.2e1 - t346;
t390 = -mrSges(7,3) / 0.2e1;
t229 = t233 ^ 2;
t389 = t229 / 0.2e1;
t227 = t231 ^ 2;
t299 = t227 + t229;
t387 = m(6) * t299;
t388 = -t387 / 0.2e1;
t386 = mrSges(5,3) - mrSges(4,2);
t385 = -Ifges(4,4) - Ifges(5,6);
t176 = mrSges(7,1) * t262 - mrSges(7,2) * t263;
t325 = t233 * mrSges(6,2);
t330 = t231 * mrSges(6,1);
t382 = t176 + t325 + t330;
t326 = t233 * mrSges(6,1);
t329 = t231 * mrSges(6,2);
t381 = t329 / 0.2e1 - t326 / 0.2e1;
t185 = -t236 * t220 + t221 * t355;
t143 = -t213 * pkin(4) + t185;
t129 = t233 * t143;
t320 = qJ(4) * t213;
t130 = t320 + t393;
t49 = -pkin(5) * t213 + t129 + (-pkin(8) * t216 - t130) * t231;
t306 = t216 * t233;
t68 = t233 * t130 + t231 * t143;
t52 = pkin(8) * t306 + t68;
t33 = t237 * t49 - t354 * t52;
t34 = t237 * t52 + t354 * t49;
t380 = t262 * t34 - t263 * t33;
t295 = t262 ^ 2 + t263 ^ 2;
t379 = -m(6) * (t231 * t59 - t233 * t60) - t231 * t168 + t233 * t170;
t255 = (-t314 / 0.2e1 + t307 / 0.2e1) * mrSges(7,3);
t249 = t255 + t346;
t332 = t254 * mrSges(7,2);
t334 = t383 * mrSges(7,1);
t265 = -t334 / 0.2e1 - t332 / 0.2e1;
t11 = t249 - t265;
t378 = t11 * qJD(1);
t309 = t263 * t148;
t313 = t262 * t151;
t360 = -t263 / 0.2e1;
t248 = (-t309 / 0.2e1 - t313 / 0.2e1) * mrSges(7,3) + t105 * t363 + t103 * t360;
t370 = mrSges(7,2) / 0.2e1;
t371 = -mrSges(7,1) / 0.2e1;
t264 = t254 * t371 + t370 * t383;
t14 = t248 + t264;
t377 = t14 * qJD(1);
t375 = -m(6) / 0.2e1;
t369 = -Ifges(7,6) / 0.2e1;
t367 = -t383 / 0.2e1;
t366 = t254 / 0.2e1;
t365 = -t173 / 0.2e1;
t208 = Ifges(7,4) * t262;
t181 = -Ifges(7,1) * t263 - t208;
t364 = t181 / 0.2e1;
t362 = -t213 / 0.2e1;
t359 = -t231 / 0.2e1;
t358 = t231 / 0.2e1;
t357 = -t233 / 0.2e1;
t356 = t233 / 0.2e1;
t175 = pkin(3) * t216 + t320;
t353 = m(5) * t175;
t352 = m(5) * t216;
t351 = mrSges(5,1) + mrSges(4,3);
t343 = Ifges(6,4) * t231;
t342 = Ifges(6,4) * t233;
t341 = Ifges(7,4) * t151;
t340 = Ifges(7,4) * t263;
t339 = Ifges(6,5) * t231;
t338 = Ifges(6,6) * t233;
t203 = t216 * mrSges(5,2);
t204 = t216 * mrSges(4,1);
t104 = mrSges(7,2) * t213 - mrSges(7,3) * t383;
t106 = -mrSges(7,1) * t213 - mrSges(7,3) * t254;
t169 = -t213 * mrSges(6,1) - t216 * t328;
t171 = t213 * mrSges(6,2) + mrSges(6,3) * t306;
t293 = t353 / 0.2e1;
t67 = -t130 * t231 + t129;
t242 = t293 + (-t231 * t67 + t233 * t68) * t374 + (-t262 * t33 - t263 * t34) * t372 + t106 * t363 + t104 * t360 + t169 * t359 + t171 * t356;
t222 = pkin(5) * t231 + qJ(4);
t243 = t293 + (-t299 * t393 - t320) * t375 + (-t172 * t383 + t173 * t254 - t222 * t213) * t373;
t282 = t227 / 0.2e1 + t389;
t331 = t216 * mrSges(6,3);
t6 = t204 - t203 + (t312 / 0.2e1 + t308 / 0.2e1) * mrSges(7,3) + t282 * t331 + (t176 / 0.2e1 + t330 / 0.2e1 + t325 / 0.2e1 + t386) * t213 + t242 + t243;
t337 = qJD(1) * t6;
t118 = -Ifges(6,6) * t213 + (Ifges(6,2) * t233 + t343) * t216;
t327 = t231 * Ifges(6,1);
t119 = -Ifges(6,5) * t213 + (t327 + t342) * t216;
t277 = -t326 + t329;
t163 = t277 * t213;
t164 = t277 * t216;
t165 = pkin(3) * t213 + t268;
t177 = -mrSges(5,2) * t213 - mrSges(5,3) * t216;
t266 = Ifges(7,5) * t366 + Ifges(7,6) * t367;
t63 = -Ifges(7,2) * t148 + t216 * Ifges(7,6) + t341;
t64 = Ifges(7,4) * t254 - Ifges(7,2) * t383 - Ifges(7,6) * t213;
t147 = Ifges(7,4) * t148;
t65 = Ifges(7,1) * t151 + t216 * Ifges(7,5) - t147;
t66 = Ifges(7,1) * t254 - Ifges(7,4) * t383 - Ifges(7,5) * t213;
t333 = t151 * mrSges(7,2);
t335 = t148 * mrSges(7,1);
t75 = t333 + t335;
t76 = t332 + t334;
t283 = -pkin(5) * t233 - pkin(4);
t98 = t216 * t283 - t184;
t99 = t213 * t283 + t185;
t1 = t284 * t204 + t98 * t75 + t99 * t76 + t34 * t103 + t30 * t104 + t33 * t105 + t29 * t106 - t148 * t64 / 0.2e1 + t63 * t367 + t151 * t66 / 0.2e1 + t65 * t366 - t141 * t163 + t143 * t164 + t67 * t168 + t59 * t169 + t68 * t170 + t60 * t171 + t175 * t177 + (-Ifges(7,5) * t151 / 0.2e1 - t148 * t369 - t284 * mrSges(4,2) + t119 * t358 + t118 * t356 + (-t339 / 0.2e1 - t338 / 0.2e1 - t385) * t213) * t213 + m(6) * (-t141 * t143 + t59 * t67 + t60 * t68) + m(7) * (t29 * t33 + t30 * t34 + t98 * t99) + ((t338 + t339 + t385) * t216 + (-Ifges(6,3) + Ifges(6,2) * t389 + Ifges(5,3) - Ifges(7,3) - Ifges(4,1) - Ifges(5,2) + Ifges(4,2) + (t342 + t327 / 0.2e1) * t231) * t213 + t266) * t216 + (t213 * mrSges(5,3) - t203 + t353) * t165;
t336 = t1 * qJD(1);
t302 = -Ifges(7,5) * t148 - Ifges(7,6) * t151;
t74 = t151 * mrSges(7,1) - mrSges(7,2) * t148;
t77 = -Ifges(7,2) * t151 - t147;
t78 = -Ifges(7,1) * t148 - t341;
t4 = t216 * t302 / 0.2e1 - t30 * t105 + t29 * t103 + t99 * t74 + (-t30 * mrSges(7,3) - t63 / 0.2e1 + t78 / 0.2e1) * t151 - (-t29 * mrSges(7,3) + t77 / 0.2e1 + t65 / 0.2e1) * t148;
t322 = t4 * qJD(1);
t5 = t254 * t103 - t383 * t105 + (t216 * t351 + t304 + t305) * t216 + (t213 * t351 - t163 - t75) * t213 + m(7) * (-t213 * t99 + t254 * t30 - t29 * t383) + m(6) * (-t143 * t213 + t216 * t275) + (m(3) * qJ(2) + mrSges(3,3)) * (t232 ^ 2 + t234 ^ 2) + (m(5) + m(4)) * (t184 * t216 - t185 * t213);
t321 = t5 * qJD(1);
t256 = t395 * t372;
t261 = m(7) * t395;
t26 = -t256 - t352 - t261 / 0.2e1 - 0.2e1 * t282 * t216 * m(6);
t319 = qJD(1) * t26;
t16 = m(7) * (t254 * t29 + t30 * t383) + t254 * t105 + t383 * t103 + (-m(5) * t165 - t177 - t379) * t216;
t316 = t16 * qJD(1);
t300 = -Ifges(7,5) * t262 + Ifges(7,6) * t263;
t174 = -mrSges(7,1) * t263 - mrSges(7,2) * t262;
t297 = t174 * qJD(6);
t296 = t375 + t373;
t289 = mrSges(7,3) * t363;
t288 = mrSges(7,3) * t361;
t278 = t387 / 0.2e1;
t273 = t231 * t68 + t233 * t67;
t102 = m(7) * t222 + mrSges(5,3) + (m(6) + m(5)) * qJ(4) + t382;
t240 = t143 * t374 + (t172 * t254 + t173 * t383 + t99) * t372 + t335 / 0.2e1 + t333 / 0.2e1 + t254 * t288 + t383 * t289 + t381 * t213;
t244 = t104 * t363 + t106 * t361 + t169 * t357 + t171 * t359 + t273 * t375 + t373 * t380;
t15 = t240 + t244;
t270 = qJD(1) * t15 + qJD(3) * t102;
t269 = qJD(1) * t74 + qJD(3) * t174;
t18 = m(7) * (-t148 * t30 - t151 * t29) - t151 * t105 - t148 * t103 + t379 * t213;
t260 = t18 * qJD(1) + qJD(4) * t294;
t253 = -t295 * t372 + t388;
t81 = t253 + t296;
t259 = qJD(1) * t294 + t81 * qJD(3);
t257 = (t309 + t313) * t372;
t42 = t257 + 0.2e1 * (t387 / 0.4e1 + m(6) / 0.4e1 + m(7) / 0.4e1) * t213;
t258 = t42 * qJD(1);
t178 = Ifges(7,2) * t263 - t208;
t179 = -Ifges(7,2) * t262 - t340;
t180 = -Ifges(7,1) * t262 + t340;
t24 = t222 * t174 - (-t179 / 0.2e1 + t180 / 0.2e1) * t263 - (t178 / 0.2e1 + t364) * t262;
t241 = -(t65 / 0.4e1 + t77 / 0.4e1) * t262 - (t78 / 0.4e1 - t63 / 0.4e1) * t263 - (t172 * t390 + t178 / 0.4e1 + t181 / 0.4e1) * t148 + (mrSges(7,3) * t365 - t179 / 0.4e1 + t180 / 0.4e1) * t151 + t172 * t103 / 0.2e1 + t105 * t365 + t216 * t300 / 0.4e1 + t222 * t74 / 0.2e1 + t99 * t174 / 0.2e1;
t247 = Ifges(7,3) * t213 / 0.2e1 + t33 * t371 + t34 * t370 - t266;
t3 = t241 + t247;
t252 = -t3 * qJD(1) - t24 * qJD(3);
t47 = m(7) * (t172 * t263 - t173 * t262) + mrSges(6,3) * t299 + mrSges(7,3) * t295 + t350 * t387;
t246 = -t141 * t374 + t216 * t381 + t98 * t372 - t265;
t8 = t255 + t246 + t392;
t250 = -t8 * qJD(1) + t47 * qJD(3);
t80 = t253 - t296;
t51 = qJD(5) * t294;
t41 = t213 * t278 + t257 + (m(6) + m(7)) * t362;
t25 = t261 / 0.2e1 + t352 / 0.2e1 - t256 + (t278 + t388 - m(5) / 0.2e1) * t216;
t13 = t248 - t264;
t12 = t249 + t265;
t10 = -t383 * t288 + t254 * t289 + t242 - t243 + t382 * t362 - t299 * t331 / 0.2e1;
t9 = -t148 * t289 + t307 * t390 + t246 - t392;
t7 = m(5) * t185 - t213 * mrSges(5,1) + t240 - t244;
t2 = t241 - t247;
t17 = [qJD(2) * t5 + qJD(3) * t1 + qJD(4) * t16 + qJD(5) * t18 + qJD(6) * t4, qJD(2) * t396 + t10 * qJD(3) + t25 * qJD(4) + t41 * qJD(5) + t12 * qJD(6) + t321, t10 * qJD(2) + t7 * qJD(4) + t9 * qJD(5) + t2 * qJD(6) + t336 + (qJ(4) * t164 + t173 * t104 + t172 * t106 + t254 * t364 + t98 * t176 + t179 * t367 + t222 * t76 + t66 * t360 + t64 * t363 + (t119 / 0.2e1 - t141 * mrSges(6,2) - t67 * mrSges(6,3) - t350 * t169) * t233 + (-t141 * mrSges(6,1) - t68 * mrSges(6,3) - t350 * t171 - t118 / 0.2e1) * t231 + (pkin(3) * mrSges(5,1) + Ifges(6,5) * t357 + Ifges(7,5) * t361 + Ifges(6,6) * t358 - t262 * t369 + Ifges(5,4) - Ifges(4,5)) * t213 + 0.2e1 * (t172 * t33 + t173 * t34 + t222 * t98) * t372 + 0.2e1 * (-qJ(4) * t141 - t273 * t350) * t374 + (-qJ(4) * mrSges(5,1) + (Ifges(6,1) * t233 - t343) * t358 + (-Ifges(6,2) * t231 + t342) * t356 + Ifges(5,5) - Ifges(4,6)) * t216 + (-m(5) * pkin(3) - mrSges(4,1) + mrSges(5,2)) * t185 + (-m(5) * qJ(4) - t386) * t184 - t380 * mrSges(7,3)) * qJD(3), t25 * qJD(2) + t7 * qJD(3) + qJD(4) * t396 + t13 * qJD(6) + t316 + t51, t41 * qJD(2) + t9 * qJD(3) + t260, t322 + t12 * qJD(2) + t2 * qJD(3) + t13 * qJD(4) + (-mrSges(7,1) * t30 - mrSges(7,2) * t29 + t302) * qJD(6); qJD(3) * t6 + qJD(4) * t26 + qJD(5) * t42 + qJD(6) * t11 - t321, 0, t337, t319, t258, -t297 + t378; -qJD(2) * t6 + qJD(4) * t15 - qJD(5) * t8 + qJD(6) * t3 - t336, -t337, qJD(4) * t102 + qJD(5) * t47 + qJD(6) * t24, qJD(5) * t80 + t270, t80 * qJD(4) + t250 (-mrSges(7,1) * t173 - mrSges(7,2) * t172 + t300) * qJD(6) - t252; -qJD(2) * t26 - qJD(3) * t15 + qJD(6) * t14 - t316 + t51, -t319, qJD(5) * t81 - t270, 0, t259, -t176 * qJD(6) + t377; -t42 * qJD(2) + t8 * qJD(3) + t74 * qJD(6) - t260, -t258, -t81 * qJD(4) - t250 + t297, -t259, 0, t269; -qJD(2) * t11 - qJD(3) * t3 - qJD(4) * t14 - qJD(5) * t74 - t322, -t378, -t174 * qJD(5) + t252, -t377, -t269, 0;];
Cq  = t17;
