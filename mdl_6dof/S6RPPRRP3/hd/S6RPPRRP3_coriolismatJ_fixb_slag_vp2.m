% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:02:46
% EndTime: 2019-03-09 02:02:53
% DurationCPUTime: 3.54s
% Computational Cost: add. (5182->392), mult. (10332->524), div. (0->0), fcn. (7919->6), ass. (0->213)
t227 = sin(qJ(5));
t225 = t227 ^ 2;
t229 = cos(qJ(5));
t226 = t229 ^ 2;
t308 = t225 + t226;
t399 = (-0.1e1 + t308) * (m(6) / 0.4e1 + m(7) / 0.4e1);
t230 = cos(qJ(4));
t356 = pkin(8) * t230;
t309 = t308 * t356;
t228 = sin(qJ(4));
t315 = t227 * t230;
t172 = -t228 * mrSges(6,2) - mrSges(6,3) * t315;
t332 = t228 * mrSges(7,3);
t176 = -mrSges(7,2) * t315 + t332;
t392 = t172 + t176;
t398 = -mrSges(6,1) - mrSges(7,1);
t397 = Ifges(7,4) + Ifges(6,5);
t396 = Ifges(7,2) + Ifges(6,3);
t395 = Ifges(6,6) - Ifges(7,6);
t207 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t283 = t207 * t227 - pkin(5);
t213 = sin(pkin(9)) * pkin(1) + qJ(3);
t355 = t228 * pkin(4);
t150 = t213 + t355 - t356;
t319 = t150 * t229;
t49 = t283 * t228 - t319;
t316 = t227 * t228;
t55 = -t207 * t316 + t319;
t394 = t49 + t55;
t393 = t230 * t228 * t399;
t312 = t229 * t230;
t174 = t228 * mrSges(6,1) - mrSges(6,3) * t312;
t175 = -t228 * mrSges(7,1) + mrSges(7,2) * t312;
t391 = t174 - t175;
t266 = pkin(5) * t229 + qJ(6) * t227;
t182 = -pkin(4) - t266;
t330 = t229 * mrSges(7,1);
t336 = t227 * mrSges(7,3);
t276 = t330 + t336;
t285 = m(7) * t182 - t276;
t348 = mrSges(7,3) * t229;
t350 = mrSges(7,1) * t227;
t275 = -t348 + t350;
t349 = mrSges(6,2) * t229;
t351 = mrSges(6,1) * t227;
t277 = t349 + t351;
t390 = t275 + t277;
t389 = -t309 / 0.2e1;
t354 = t228 * pkin(8);
t190 = pkin(4) * t230 + t354;
t313 = t229 * t190;
t65 = -t207 * t315 + t313;
t66 = t227 * t190 + t207 * t312;
t260 = -t65 * t227 + t66 * t229;
t57 = qJ(6) * t230 + t66;
t58 = t283 * t230 - t313;
t261 = t227 * t58 + t229 * t57;
t296 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t383 = t230 ^ 2;
t384 = t228 ^ 2;
t388 = t383 / 0.2e1 + t384 / 0.2e1;
t223 = m(7) * qJ(6) + mrSges(7,3);
t338 = t226 * mrSges(6,3);
t339 = t226 * mrSges(7,2);
t340 = t225 * mrSges(6,3);
t341 = t225 * mrSges(7,2);
t352 = m(7) * qJD(6);
t387 = (-mrSges(5,2) + t338 + t339 + t340 + t341) * qJD(4) + t227 * t352;
t364 = t229 / 0.2e1;
t365 = -t229 / 0.2e1;
t371 = -t227 / 0.2e1;
t386 = t174 * t365 + t175 * t364 + t371 * t392;
t378 = m(7) / 0.2e1;
t382 = -m(6) / 0.2e1;
t318 = t207 * t229;
t320 = t150 * t227;
t48 = t320 + (qJ(6) + t318) * t228;
t314 = t228 * t229;
t56 = t207 * t314 + t320;
t265 = pkin(5) * t227 - qJ(6) * t229;
t254 = -t207 + t265;
t67 = t254 * t228;
t385 = (t227 * t49 + t229 * t48 + t67) * t378 + (t227 * t55 - t229 * t56) * t382;
t381 = m(6) / 0.2e1;
t379 = -m(7) / 0.2e1;
t376 = mrSges(6,1) / 0.2e1;
t375 = -mrSges(7,1) / 0.2e1;
t374 = qJ(6) / 0.2e1;
t331 = t229 * mrSges(6,1);
t337 = t227 * mrSges(6,2);
t278 = t331 - t337;
t153 = t278 * t230;
t373 = -t153 / 0.2e1;
t154 = t230 * t275;
t372 = -t154 / 0.2e1;
t369 = t227 / 0.2e1;
t368 = t227 / 0.4e1;
t367 = t228 / 0.2e1;
t362 = -t230 / 0.2e1;
t361 = t230 / 0.2e1;
t360 = m(6) * t207;
t151 = t266 * t230;
t359 = m(7) * t151;
t357 = m(7) * t229;
t353 = Ifges(7,5) - Ifges(6,4);
t347 = Ifges(6,4) * t227;
t346 = Ifges(6,4) * t229;
t345 = Ifges(7,5) * t227;
t344 = Ifges(7,5) * t229;
t343 = Ifges(6,6) * t228;
t342 = pkin(8) * qJD(4);
t335 = t227 * t48;
t334 = t227 * t56;
t328 = t230 * mrSges(7,1);
t325 = -t278 - mrSges(5,1);
t263 = t55 * t229 + t334;
t264 = -t49 * t229 + t335;
t12 = t228 * mrSges(5,1) + t230 * mrSges(5,2) + mrSges(4,3) + t391 * t229 + t392 * t227 + m(7) * t264 + m(6) * t263 + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t213;
t324 = qJD(1) * t12;
t68 = t254 * t230;
t23 = m(7) * (t48 * t228 - t68 * t312) - t154 * t312 + t228 * t176;
t323 = qJD(1) * t23;
t233 = (t263 - t264) * t378 + t386;
t152 = t276 * t230;
t253 = -t359 / 0.2e1 - t152 / 0.2e1 + t373;
t232 = t253 * t230 + (t233 + (mrSges(7,2) + mrSges(6,3)) * t230 * (-t225 / 0.2e1 - t226 / 0.2e1)) * t228;
t235 = t266 * t379 + t337 / 0.2e1 - t336 / 0.2e1 - t331 / 0.2e1 - t330 / 0.2e1;
t10 = t232 + t235;
t322 = t10 * qJD(1);
t11 = -t253 * t228 + ((-t341 / 0.2e1 - t339 / 0.2e1 - t338 / 0.2e1 - t340 / 0.2e1) * t230 + t233) * t230;
t321 = t11 * qJD(1);
t317 = t207 * t230;
t74 = (0.1e1 / 0.2e1 + t388) * t357;
t310 = t74 * qJD(1);
t307 = qJD(4) * t228;
t306 = qJD(4) * t230;
t305 = qJD(5) * t227;
t304 = qJD(5) * t229;
t303 = qJD(5) * t230;
t300 = t58 * t378;
t297 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t289 = -t312 / 0.2e1;
t288 = -t174 / 0.2e1 + t175 / 0.2e1;
t287 = t176 / 0.2e1 + t172 / 0.2e1;
t284 = -mrSges(7,2) * qJ(6) - Ifges(6,6);
t280 = -m(6) * pkin(4) + t325;
t279 = mrSges(7,2) * pkin(5) - t397;
t274 = t229 * Ifges(6,1) - t347;
t273 = Ifges(6,1) * t227 + t346;
t272 = t229 * Ifges(7,1) + t345;
t271 = Ifges(7,1) * t227 - t344;
t270 = -Ifges(6,2) * t227 + t346;
t269 = Ifges(6,2) * t229 + t347;
t268 = Ifges(7,3) * t227 + t344;
t267 = -Ifges(7,3) * t229 + t345;
t209 = Ifges(7,6) * t312;
t210 = Ifges(7,5) * t312;
t249 = t272 * t230;
t128 = Ifges(7,4) * t228 + t249;
t250 = t274 * t230;
t130 = Ifges(6,5) * t228 + t250;
t241 = t55 * mrSges(6,3) - t297 * t228 - t128 / 0.2e1 - t130 / 0.2e1;
t124 = Ifges(7,6) * t228 + Ifges(7,3) * t315 + t210;
t248 = t270 * t230;
t126 = t248 + t343;
t256 = t56 * mrSges(6,3) - t124 / 0.2e1 + t126 / 0.2e1;
t5 = t68 * t152 + t209 * t367 + (-t207 * t153 + (-t48 * mrSges(7,2) - t343 / 0.2e1 + t210 / 0.2e1 + Ifges(6,4) * t289 - t256) * t229 + (-t49 * mrSges(7,2) + (-Ifges(7,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t315 + (Ifges(7,3) / 0.2e1 - Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(7,1) / 0.2e1) * t312 + t241) * t227) * t230 + (m(7) * t68 + t154) * t151 + (m(7) * t49 - t391) * t56 + (m(7) * t48 + t392) * t55;
t259 = t5 * qJD(1) + t11 * qJD(2);
t208 = mrSges(7,2) * t314;
t238 = (-t227 * t68 + (-t182 * t230 + t354) * t229) * t379 + t154 * t369;
t21 = -t208 + (-t276 * t364 + t375) * t230 + t300 + t238;
t61 = t285 * t227;
t258 = qJD(1) * t21 + qJD(4) * t61;
t29 = t332 + 0.2e1 * (t320 / 0.4e1 - t56 / 0.4e1 + (t318 / 0.4e1 + t374) * t228) * m(7);
t257 = qJD(1) * t29 + qJD(5) * t223;
t252 = m(6) * t260;
t251 = m(7) * t265;
t123 = t230 * Ifges(7,6) - t268 * t228;
t125 = t230 * Ifges(6,6) - t270 * t228;
t127 = t230 * Ifges(7,4) - t272 * t228;
t129 = t230 * Ifges(6,5) - t274 * t228;
t173 = -t208 - t328;
t177 = mrSges(7,2) * t316 + t230 * mrSges(7,3);
t3 = -t67 * t154 + t66 * t172 + t49 * t173 + t65 * t174 + t58 * t175 + t57 * t176 + t48 * t177 + m(7) * (t48 * t57 + t49 * t58 - t67 * t68) + m(6) * (t55 * t65 + t56 * t66) + (t213 * mrSges(5,1) + t55 * mrSges(6,1) - t56 * mrSges(6,2) - Ifges(5,4) * t230 + (t127 / 0.2e1 + t129 / 0.2e1 + t297 * t230) * t229 + (t123 / 0.2e1 - t125 / 0.2e1 + t296 * t230) * t227) * t230 + (-t213 * mrSges(5,2) + Ifges(5,4) * t228 + (t68 * mrSges(7,3) + t241) * t229 + (-t68 * mrSges(7,1) - t296 * t228 + t256) * t227 + (-Ifges(5,1) + Ifges(5,2) + (0.2e1 * t349 + 0.2e1 * t351 - t360) * t207 + t396) * t230) * t228;
t237 = t372 + (t261 + t68) * t379 + t173 * t371 + t177 * t365;
t6 = (((-mrSges(7,3) / 0.2e1 + mrSges(6,2) / 0.2e1) * t228 + t287) * t229 + ((mrSges(7,1) / 0.2e1 + t376) * t228 + t288) * t227 - t228 * t360 / 0.2e1 + t385) * t228 + (-t252 / 0.2e1 + t317 * t381 + t237) * t230;
t8 = (t288 * t227 + t287 * t229 + t385) * t230 + ((t260 - 0.2e1 * t317) * t381 - t237 + t390 * t361) * t228;
t246 = t3 * qJD(1) - t6 * qJD(2) + t8 * qJD(3);
t245 = t277 * t362;
t34 = 0.2e1 * (t383 - t384) * t399;
t37 = 0.4e1 * t393;
t244 = t8 * qJD(1) + t34 * qJD(2) + t37 * qJD(3);
t38 = -0.4e1 * t393;
t243 = -t6 * qJD(1) + t38 * qJD(2) + t34 * qJD(3);
t231 = t269 * t289 + t207 * t245 + (t182 * t151 + t265 * t68) * t378 + pkin(4) * t373 + t182 * t152 / 0.2e1 - t151 * t276 / 0.2e1 + t68 * t275 / 0.2e1 + t265 * t154 / 0.2e1 + t267 * t312 / 0.2e1 - (t248 + t126) * t227 / 0.4e1 + (t210 + t124) * t368 + (-t227 * t395 + t229 * t397) * t228 / 0.4e1 + t389 * mrSges(6,3) + (-Ifges(7,1) * t368 - t273 / 0.2e1 - t271 / 0.4e1 + t268 / 0.4e1) * t315 + (t128 + t250 + t249 + t130) * t229 / 0.4e1 + ((t394 * t229 + (-t48 + t56) * t227) * t378 + t386) * pkin(8) + (t394 * t364 + t334 / 0.2e1 - t335 / 0.2e1 + t389) * mrSges(7,2);
t234 = (-pkin(5) * t58 + qJ(6) * t57) * t378 - pkin(5) * t173 / 0.2e1 + t177 * t374 + t57 * mrSges(7,3) / 0.2e1 + t58 * t375 + t65 * t376 - t66 * mrSges(6,2) / 0.2e1;
t1 = -t231 + t234 + t396 * t361 - t296 * t316 - t397 * t314 / 0.2e1;
t18 = t265 * t276 + (pkin(4) * mrSges(6,2) + t353 * t229) * t229 + (pkin(4) * mrSges(6,1) - t353 * t227 + (-Ifges(6,1) - Ifges(7,1) + Ifges(6,2) + Ifges(7,3)) * t229) * t227 + (-t251 - t275) * t182;
t242 = t1 * qJD(1) + t18 * qJD(4);
t239 = -t351 / 0.2e1 - t350 / 0.2e1 - t349 / 0.2e1 + t348 / 0.2e1;
t236 = (-m(7) * t266 - t276 - t278) * qJD(5);
t183 = (m(7) * pkin(8) + mrSges(7,2)) * t229;
t75 = (-0.1e1 / 0.2e1 + t388) * t357;
t33 = t34 * qJD(4);
t28 = t176 + (t320 + (0.2e1 * qJ(6) + t318) * t228 + t56) * t378;
t27 = (t251 - t239) * t228 + t390 * t367;
t26 = t251 * t362 + t245 + t372 + (t265 * t379 + t239) * t230;
t22 = -t276 * t289 + t300 - t328 / 0.2e1 - t238;
t9 = t232 - t235;
t7 = t8 * qJD(4);
t4 = -qJD(4) * t6 + qJD(5) * t11;
t2 = (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t230 + t231 + (-t296 * t227 - t297 * t229) * t228 + t234;
t13 = [qJD(3) * t12 + qJD(4) * t3 + qJD(5) * t5 + qJD(6) * t23, t4, qJD(5) * t9 + qJD(6) * t75 + t324 + t7 (-mrSges(5,2) * t317 - Ifges(5,6) * t230 + t123 * t365 + t125 * t364 - t285 * t67 + (t127 + t129) * t369 + (t227 * t397 + t229 * t395) * t361 + t260 * mrSges(6,3) + t261 * mrSges(7,2)) * qJD(4) + t2 * qJD(5) + t22 * qJD(6) + ((-t230 * mrSges(6,2) + t177) * t229 + (-t230 * mrSges(6,1) + t173) * t227 + m(7) * t261 + t252) * t342 + (pkin(4) * t277 - t182 * t275 + t280 * t207 + t267 * t371 + t269 * t369 - Ifges(5,5) + (t271 + t273) * t365) * t307 + t246, t9 * qJD(3) + t2 * qJD(4) + (t209 + (-m(7) * pkin(5) + t398) * t56 + (-mrSges(6,2) + t223) * t55) * qJD(5) + t28 * qJD(6) + (t227 * t279 + t229 * t284) * t303 + t259, qJD(3) * t75 + qJD(4) * t22 + qJD(5) * t28 + t323; t4, t38 * qJD(4), t33, t27 * qJD(5) + (t280 + t285) * t306 + (0.2e1 * (t379 + t382) * t342 * t308 - t387) * t228 + t243, t321 + t27 * qJD(4) - qJD(5) * t359 + ((mrSges(6,2) - mrSges(7,3)) * t305 + (qJD(5) * t398 + t352) * t229) * t230 (-t227 * t307 + t229 * t303) * m(7); qJD(5) * t10 + qJD(6) * t74 - t324 + t7, t33, t37 * qJD(4), t26 * qJD(5) + (-t276 + t325) * t307 + 0.2e1 * ((t309 - t355) * t381 + (t182 * t228 + t309) * t378) * qJD(4) + t387 * t230 + t244, t322 + t26 * qJD(4) + (t229 * t352 + t236) * t228, t310 + (t227 * t306 + t228 * t304) * m(7); -qJD(5) * t1 - qJD(6) * t21 - t246, -t243, -t244, -qJD(5) * t18 - qJD(6) * t61, t183 * qJD(6) - t279 * t304 + (Ifges(7,6) + t284) * t305 + pkin(8) * t236 - t242, qJD(5) * t183 - t258; -qJD(3) * t10 + qJD(4) * t1 + qJD(6) * t29 - t259, -t321, -t322, t242, t223 * qJD(6), t257; -qJD(3) * t74 + qJD(4) * t21 - qJD(5) * t29 - t323, 0, -t310, t258, -t257, 0;];
Cq  = t13;
