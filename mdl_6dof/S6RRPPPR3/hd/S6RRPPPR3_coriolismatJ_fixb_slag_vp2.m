% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPPPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:55
% EndTime: 2019-03-09 08:14:01
% DurationCPUTime: 3.65s
% Computational Cost: add. (8604->419), mult. (15316->550), div. (0->0), fcn. (14566->6), ass. (0->217)
t250 = sin(pkin(9));
t251 = cos(pkin(9));
t361 = sin(qJ(6));
t362 = cos(qJ(6));
t201 = t362 * t250 + t361 * t251;
t137 = t201 ^ 2;
t200 = t361 * t250 - t362 * t251;
t138 = t200 ^ 2;
t313 = -t138 - t137;
t391 = m(7) * (t137 / 0.2e1 + t313 / 0.2e1 + t138 / 0.2e1);
t254 = cos(qJ(2));
t171 = t200 * t254;
t322 = t171 * t201;
t169 = t201 * t254;
t328 = t169 * t200;
t380 = m(7) / 0.2e1;
t306 = (t322 - t328) * t380;
t392 = mrSges(7,1) / 0.2e1;
t372 = -t200 / 0.2e1;
t369 = -t201 / 0.2e1;
t390 = t169 * mrSges(7,3);
t194 = Ifges(7,4) * t200;
t133 = -Ifges(7,1) * t201 + t194;
t389 = Ifges(7,2) * t201 + t133 + t194;
t248 = t250 ^ 2;
t249 = t251 ^ 2;
t309 = -t249 - t248;
t240 = t254 * qJ(3);
t255 = -pkin(2) - pkin(3);
t247 = -qJ(5) + t255;
t253 = sin(qJ(2));
t149 = pkin(4) * t254 + t247 * t253 + t240;
t239 = t254 * qJ(4);
t218 = t254 * pkin(7) - t239;
t111 = t251 * t149 - t218 * t250;
t315 = t251 * t253;
t80 = pkin(5) * t254 - pkin(8) * t315 + t111;
t112 = t250 * t149 + t251 * t218;
t318 = t250 * t253;
t85 = -pkin(8) * t318 + t112;
t44 = -t361 * t85 + t362 * t80;
t45 = t361 * t80 + t362 * t85;
t388 = t200 * t45 + t201 * t44;
t387 = Ifges(3,4) + Ifges(5,4) - Ifges(4,5);
t145 = -mrSges(7,2) * t253 + t390;
t147 = mrSges(7,1) * t253 - t171 * mrSges(7,3);
t323 = t171 * t200;
t327 = t169 * t201;
t371 = t200 / 0.2e1;
t266 = (t323 / 0.2e1 + t327 / 0.2e1) * mrSges(7,3) + t147 * t371 + t145 * t369;
t168 = t201 * t253;
t170 = t200 * t253;
t379 = -mrSges(7,2) / 0.2e1;
t278 = -t168 * t379 + t170 * t392;
t11 = t266 + t278;
t386 = -t11 * qJD(1) - qJD(4) * t391;
t271 = (-t322 / 0.2e1 + t328 / 0.2e1) * mrSges(7,3);
t276 = t145 * t372 + t147 * t369;
t264 = t271 + t276;
t347 = t170 * mrSges(7,2);
t349 = t168 * mrSges(7,1);
t279 = -t349 / 0.2e1 + t347 / 0.2e1;
t13 = t264 + t279;
t385 = -t13 * qJD(1) - qJD(3) * t391;
t199 = t255 * t253 + t240;
t384 = m(5) / 0.2e1;
t383 = -m(6) / 0.2e1;
t382 = m(6) / 0.2e1;
t381 = -m(7) / 0.2e1;
t377 = -t168 / 0.2e1;
t376 = t169 / 0.2e1;
t374 = -t170 / 0.2e1;
t373 = t171 / 0.2e1;
t368 = t201 / 0.2e1;
t367 = -t250 / 0.2e1;
t366 = t250 / 0.2e1;
t365 = -t251 / 0.2e1;
t364 = t251 / 0.2e1;
t360 = m(5) * t199;
t359 = -mrSges(4,2) + mrSges(5,3);
t358 = pkin(8) - t247;
t252 = qJ(3) + pkin(4);
t355 = Ifges(6,4) * t250;
t354 = Ifges(6,4) * t251;
t353 = Ifges(7,4) * t171;
t352 = Ifges(7,4) * t201;
t207 = t358 * t250;
t208 = t358 * t251;
t126 = t207 * t362 + t208 * t361;
t127 = t207 * t361 - t208 * t362;
t131 = -t200 * mrSges(7,1) - t201 * mrSges(7,2);
t214 = t251 * mrSges(6,1) - t250 * mrSges(6,2);
t230 = pkin(5) * t251 + t252;
t277 = -t168 * t369 - t170 * t372;
t298 = t248 / 0.2e1 + t249 / 0.2e1;
t292 = t298 * t253;
t295 = t309 * t247;
t259 = (-t131 / 0.2e1 - t214 / 0.2e1) * t254 + t277 * mrSges(7,3) + mrSges(6,3) * t292 - t199 * t384 + (-t252 * t254 + t253 * t295) * t382 + (t126 * t168 + t127 * t170 - t230 * t254) * t380;
t144 = -mrSges(7,2) * t254 - t168 * mrSges(7,3);
t146 = mrSges(7,1) * t254 + t170 * mrSges(7,3);
t335 = t254 * mrSges(6,2);
t209 = -mrSges(6,3) * t318 - t335;
t336 = t254 * mrSges(6,1);
t211 = -mrSges(6,3) * t315 + t336;
t261 = t360 / 0.2e1 + (t251 * t111 + t250 * t112) * t382 + (-t200 * t44 + t201 * t45) * t380 + t146 * t372 + t144 * t368 + t209 * t366 + t211 * t364;
t310 = t254 * mrSges(5,1) + t253 * mrSges(5,2);
t5 = t259 - t261 - t310;
t351 = qJD(1) * t5;
t213 = -t254 * pkin(2) - t253 * qJ(3) - pkin(1);
t196 = t254 * pkin(3) - t213;
t148 = pkin(4) * t253 + qJ(5) * t254 + t196;
t237 = t253 * qJ(4);
t215 = pkin(7) * t253 - t237;
t100 = t250 * t148 + t251 * t215;
t342 = t250 * Ifges(6,2);
t166 = Ifges(6,6) * t254 + (-t342 + t354) * t253;
t167 = Ifges(6,5) * t254 + (t251 * Ifges(6,1) - t355) * t253;
t303 = pkin(5) * t250 - pkin(7);
t182 = t253 * t303 + t237;
t183 = -t254 * t303 - t239;
t340 = t251 * mrSges(6,2);
t343 = t250 * mrSges(6,1);
t290 = t340 + t343;
t187 = t290 * t253;
t188 = t290 * t254;
t317 = t250 * t254;
t337 = t253 * mrSges(6,2);
t210 = mrSges(6,3) * t317 - t337;
t314 = t251 * t254;
t338 = t253 * mrSges(6,1);
t212 = mrSges(6,3) * t314 + t338;
t219 = mrSges(5,1) * t253 - mrSges(5,2) * t254;
t280 = Ifges(7,5) * t374 + Ifges(7,6) * t377;
t297 = m(4) * t213 - t254 * mrSges(4,1) - t253 * mrSges(4,3);
t339 = t251 * Ifges(6,5);
t341 = t250 * Ifges(6,6);
t99 = t251 * t148 - t215 * t250;
t75 = pkin(5) * t253 + pkin(8) * t314 + t99;
t83 = pkin(8) * t317 + t100;
t42 = -t361 * t83 + t362 * t75;
t43 = t361 * t75 + t362 * t83;
t90 = -Ifges(7,4) * t170 - Ifges(7,2) * t168 + Ifges(7,6) * t254;
t91 = Ifges(7,2) * t169 + t253 * Ifges(7,6) + t353;
t92 = -Ifges(7,1) * t170 - Ifges(7,4) * t168 + Ifges(7,5) * t254;
t165 = Ifges(7,4) * t169;
t93 = Ifges(7,1) * t171 + t253 * Ifges(7,5) + t165;
t95 = -t347 + t349;
t346 = t171 * mrSges(7,2);
t348 = t169 * mrSges(7,1);
t96 = t346 - t348;
t1 = t43 * t144 + t45 * t145 + t42 * t146 + t44 * t147 + t91 * t377 + t90 * t376 + t93 * t374 + t92 * t373 + t182 * t96 + t183 * t95 + t100 * t209 + t112 * t210 + t99 * t211 + t111 * t212 + t215 * t188 + t218 * t187 + t199 * t219 + t297 * (pkin(2) * t253 - t240) + m(6) * (t100 * t112 + t111 * t99 - t215 * t218) + m(7) * (t182 * t183 + t42 * t44 + t43 * t45) + (-pkin(1) * mrSges(3,2) - t213 * mrSges(4,3) + Ifges(7,5) * t373 + Ifges(7,6) * t376 + t166 * t366 + t167 * t365 + (-t339 / 0.2e1 + t341 / 0.2e1 + t387) * t254) * t254 + (-pkin(1) * mrSges(3,1) + t213 * mrSges(4,1) + (-t249 * Ifges(6,1) / 0.2e1 - Ifges(3,2) + Ifges(6,3) - Ifges(5,1) + Ifges(5,2) + Ifges(4,1) - Ifges(4,3) + Ifges(3,1) + Ifges(7,3) + (t354 - t342 / 0.2e1) * t250) * t254 + t280 + (t339 - t341 - t387) * t253) * t253 + (t310 + t360) * t196;
t350 = t1 * qJD(1);
t312 = Ifges(7,5) * t169 - Ifges(7,6) * t171;
t94 = t171 * mrSges(7,1) + mrSges(7,2) * t169;
t97 = -Ifges(7,2) * t171 + t165;
t98 = Ifges(7,1) * t169 - t353;
t4 = t42 * t145 - t43 * t147 + t253 * t312 / 0.2e1 + t183 * t94 - (t43 * mrSges(7,3) - t98 / 0.2e1 + t91 / 0.2e1) * t171 + (-t42 * mrSges(7,3) + t93 / 0.2e1 + t97 / 0.2e1) * t169;
t334 = t4 * qJD(1);
t31 = m(5) * t253 + 0.2e1 * m(6) * t292 + 0.2e1 * m(7) * t277;
t333 = qJD(1) * t31;
t287 = -t251 * t100 + t250 * t99;
t316 = t251 * t210;
t319 = t250 * t212;
t320 = t218 * t254;
t15 = t170 * t145 + t168 * t147 + (t254 * mrSges(5,3) + t188 - t96) * t254 + (t253 * mrSges(5,3) - t316 + t319) * t253 + m(7) * (t168 * t42 + t170 * t43 - t183 * t254) + m(6) * (t253 * t287 - t320) + m(5) * (-t215 * t253 - t320);
t330 = t15 * qJD(1);
t267 = m(6) * (t100 * t250 + t251 * t99) + t251 * t212 + t250 * t210;
t17 = m(7) * (t168 * t43 - t170 * t42) + t168 * t145 - t170 * t147 + (m(5) * t196 + t219 + t267 - t297) * t253;
t326 = t17 * qJD(1);
t324 = t170 * t201;
t128 = t200 * t168;
t311 = Ifges(7,5) * t200 + Ifges(7,6) * t201;
t130 = -t201 * mrSges(7,1) + mrSges(7,2) * t200;
t308 = t130 * qJD(6);
t307 = t383 + t381;
t294 = t309 * t383;
t293 = t309 * t382;
t291 = m(7) * (-t128 + t324);
t289 = Ifges(7,1) * t200 + t352;
t132 = Ifges(7,2) * t200 - t352;
t258 = -t127 * t147 / 0.2e1 + t183 * t130 / 0.2e1 + t230 * t94 / 0.2e1 + t253 * t311 / 0.4e1 + t389 * t169 / 0.4e1 + (t93 + t97) * t200 / 0.4e1 + (t91 / 0.4e1 - t98 / 0.4e1) * t201 + (t145 / 0.2e1 - t390 / 0.2e1) * t126 - (t127 * mrSges(7,3) / 0.2e1 + t132 / 0.4e1 - t289 / 0.4e1) * t171;
t265 = Ifges(7,3) * t254 / 0.2e1 + t44 * t392 + t45 * t379 + t280;
t2 = t258 - t265;
t21 = -t230 * t130 + t132 * t369 + t289 * t368 + t372 * t389;
t288 = t2 * qJD(1) - t21 * qJD(2);
t262 = (t128 / 0.2e1 - t324 / 0.2e1) * mrSges(7,3) + (-t126 * t170 + t127 * t168 + t183) * t380 - t348 / 0.2e1 + t346 / 0.2e1;
t285 = -t250 * t111 + t251 * t112;
t263 = t144 * t371 + t146 * t368 + t285 * t383 - t381 * t388;
t10 = (-t335 / 0.2e1 - t209 / 0.2e1) * t251 + (-t336 / 0.2e1 + t211 / 0.2e1) * t250 + t218 * t382 + t262 + t263;
t60 = mrSges(5,1) + mrSges(4,3) + (m(5) + m(4)) * qJ(3) + m(6) * t252 + m(7) * t230 + t214 + t131;
t286 = qJD(1) * t10 + qJD(2) * t60;
t283 = qJD(1) * t94 + qJD(2) * t130;
t20 = m(7) * (t169 * t43 - t171 * t42) + t169 * t145 - t171 * t147 + t267 * t254;
t275 = t20 * qJD(1) + qJD(3) * t306;
t272 = (t323 + t327) * t380;
t47 = t272 + (t380 + (0.1e1 / 0.2e1 + t298) * m(6)) * t254;
t274 = t47 * qJD(1);
t270 = t313 * t380 + t293;
t55 = t270 + t307;
t273 = qJD(1) * t306 + t55 * qJD(2);
t33 = -m(7) * (t126 * t201 + t127 * t200) - m(6) * t295 + t313 * mrSges(7,3) + t309 * mrSges(6,3);
t260 = t271 + t287 * t382 + (-t126 * t171 + t127 * t169 + t200 * t43 + t201 * t42) * t380 - t276;
t268 = t182 * t381 + t215 * t382 + t279;
t9 = (-t210 / 0.2e1 - t337 / 0.2e1) * t251 + (t212 / 0.2e1 - t338 / 0.2e1) * t250 + t260 + t268;
t269 = t9 * qJD(1) - t33 * qJD(2);
t56 = qJD(5) * t306;
t54 = t270 - t307;
t52 = qJD(6) * t391;
t48 = t254 * t294 + t272 - (m(6) + m(7)) * t254 / 0.2e1;
t32 = (t293 + t294) * t253;
t14 = t264 - t279;
t12 = t266 - t278;
t8 = t319 / 0.2e1 - t316 / 0.2e1 + mrSges(6,2) * t315 / 0.2e1 + mrSges(6,1) * t318 / 0.2e1 + t260 - t268;
t7 = t211 * t367 + t209 * t364 + (-t340 / 0.2e1 - t343 / 0.2e1 + m(4) * pkin(7) - t359) * t254 + 0.2e1 * (m(6) / 0.4e1 + t384) * t218 + t262 - t263;
t6 = t259 + t261;
t3 = t258 + t265;
t16 = [qJD(2) * t1 + qJD(3) * t17 + qJD(4) * t15 + qJD(5) * t20 + qJD(6) * t4, t7 * qJD(3) + t6 * qJD(4) + t8 * qJD(5) + t3 * qJD(6) + t350 + (-t215 * mrSges(5,1) + t218 * mrSges(5,2) + t126 * t146 + t127 * t144 + t182 * t131 + t132 * t377 + t133 * t374 + t252 * t187 - t215 * t214 + t230 * t95 + t92 * t369 + t90 * t371 + (t247 * t209 - t112 * mrSges(6,3) - t166 / 0.2e1) * t251 + (-t247 * t211 + t111 * mrSges(6,3) - t167 / 0.2e1) * t250 + 0.2e1 * (-qJ(3) * t215 + t218 * t255) * t384 + 0.2e1 * (-t215 * t252 + t247 * t285) * t382 + 0.2e1 * (t126 * t44 + t127 * t45 + t182 * t230) * t380 + (-pkin(2) * mrSges(4,2) - t255 * mrSges(5,3) + Ifges(5,6) + Ifges(4,4) + Ifges(3,5) + Ifges(7,5) * t369 + Ifges(7,6) * t371 + Ifges(6,5) * t367 + Ifges(6,6) * t365 + (-m(4) * pkin(2) - mrSges(3,1) - mrSges(4,1)) * pkin(7)) * t254 + ((-Ifges(6,1) * t250 - t354) * t364 + (-Ifges(6,2) * t251 - t355) * t367 + Ifges(4,6) - Ifges(3,6) - Ifges(5,5) + t359 * qJ(3) + (-m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3)) * pkin(7)) * t253 + t388 * mrSges(7,3)) * qJD(2), t7 * qJD(2) + qJD(3) * t291 + t32 * qJD(4) + t12 * qJD(6) + t326 + t56, t6 * qJD(2) + t32 * qJD(3) + qJD(4) * t291 + t48 * qJD(5) + t14 * qJD(6) + t330, t8 * qJD(2) + t48 * qJD(4) + t275, t334 + t3 * qJD(2) + t12 * qJD(3) + t14 * qJD(4) + (-mrSges(7,1) * t43 - mrSges(7,2) * t42 + t312) * qJD(6); qJD(3) * t10 + qJD(4) * t5 + qJD(5) * t9 + qJD(6) * t2 - t350, qJD(3) * t60 - qJD(5) * t33 - qJD(6) * t21, qJD(5) * t54 + t286, t351, t54 * qJD(3) + t269 (-mrSges(7,1) * t127 - mrSges(7,2) * t126 + t311) * qJD(6) + t288; -qJD(2) * t10 - qJD(4) * t31 + qJD(6) * t11 - t326 + t56, qJD(5) * t55 - t286, 0, t52 - t333, t273, -t131 * qJD(6) - t386; -qJD(2) * t5 + qJD(3) * t31 + qJD(5) * t47 + qJD(6) * t13 - t330, -t351, t52 + t333, 0, t274, t308 - t385; -t9 * qJD(2) - t47 * qJD(4) + t94 * qJD(6) - t275, -t55 * qJD(3) - t269 + t308, -t273, -t274, 0, t283; -qJD(2) * t2 - qJD(3) * t11 - qJD(4) * t13 - qJD(5) * t94 - t334, -qJD(5) * t130 - t288, t386, t385, -t283, 0;];
Cq  = t16;
