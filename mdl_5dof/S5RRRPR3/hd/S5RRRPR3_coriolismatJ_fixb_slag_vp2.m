% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR3
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:54
% EndTime: 2020-01-03 12:09:02
% DurationCPUTime: 3.78s
% Computational Cost: add. (11144->249), mult. (21423->324), div. (0->0), fcn. (22873->8), ass. (0->160)
t218 = sin(qJ(3));
t219 = sin(qJ(2));
t328 = t219 * pkin(1);
t210 = pkin(7) + t328;
t277 = qJ(4) + t210;
t182 = t277 * t218;
t221 = cos(qJ(3));
t183 = t277 * t221;
t215 = sin(pkin(9));
t216 = cos(pkin(9));
t142 = -t215 * t182 + t216 * t183;
t189 = -t215 * t218 + t216 * t221;
t180 = t189 * pkin(8);
t104 = t180 + t142;
t217 = sin(qJ(5));
t220 = cos(qJ(5));
t190 = -t215 * t221 - t216 * t218;
t181 = t190 * pkin(8);
t346 = -t216 * t182 - t215 * t183;
t352 = t346 + t181;
t374 = -t104 * t217 + t220 * t352;
t66 = t104 * t220 + t217 * t352;
t390 = -t66 * mrSges(6,1) - t374 * mrSges(6,2);
t319 = -qJ(4) - pkin(7);
t198 = t319 * t218;
t200 = t319 * t221;
t159 = t215 * t198 - t216 * t200;
t125 = t180 + t159;
t345 = t216 * t198 + t215 * t200;
t353 = t345 + t181;
t373 = -t125 * t217 + t220 * t353;
t89 = t125 * t220 + t217 * t353;
t389 = -t89 * mrSges(6,1) - t373 * mrSges(6,2);
t153 = t189 * t217 - t190 * t220;
t257 = t220 * t189 + t190 * t217;
t378 = Ifges(6,5) * t257 - Ifges(6,6) * t153;
t8 = t378 + t390;
t388 = t8 * qJD(5);
t11 = t378 + t389;
t387 = t11 * qJD(5);
t214 = t221 ^ 2;
t380 = t221 * (Ifges(4,1) - Ifges(4,2));
t375 = Ifges(5,1) - Ifges(5,2);
t368 = t189 / 0.2e1;
t372 = Ifges(6,1) - Ifges(6,2);
t261 = -t221 * pkin(3) - pkin(2);
t164 = -t189 * pkin(4) + t261;
t222 = cos(qJ(2));
t332 = pkin(1) * t222;
t160 = t164 - t332;
t350 = t257 * mrSges(6,2);
t355 = t153 * mrSges(6,1);
t358 = t355 + t350;
t252 = t160 * t358;
t84 = t164 * t358;
t148 = Ifges(6,4) * t257;
t357 = Ifges(6,4) * t153;
t255 = Ifges(6,1) * t257 - t357;
t154 = -mrSges(5,1) * t189 - mrSges(5,2) * t190;
t329 = pkin(3) * t218;
t256 = -pkin(4) * t190 + t329;
t91 = -mrSges(6,1) * t257 + mrSges(6,2) * t153;
t290 = t154 * t329 + t256 * t91;
t308 = Ifges(6,2) * t257;
t309 = Ifges(5,4) * t190;
t310 = Ifges(5,4) * t189;
t366 = t257 / 0.2e1;
t371 = t290 + 0.2e1 * t366 * t148 + 0.2e1 * t310 * t368 + Ifges(4,4) * t214 + (-t308 / 0.2e1 - t357 / 0.2e1 + t255 / 0.2e1 + t372 * t366) * t153 + (-Ifges(4,4) * t218 + t380) * t218 + (-0.2e1 * t375 * t368 - t309) * t190;
t267 = t355 / 0.2e1;
t365 = m(5) * t329;
t213 = t218 ^ 2;
t360 = t213 + t214;
t259 = -t190 * mrSges(5,1) + t189 * mrSges(5,2);
t130 = t261 * t259;
t197 = t261 - t332;
t246 = t164 * t256;
t247 = t160 * t256;
t251 = t197 * t259;
t295 = t221 * mrSges(4,2);
t296 = t218 * mrSges(4,1);
t201 = t295 + t296;
t211 = -pkin(2) - t332;
t283 = t211 * t201;
t301 = t257 * mrSges(6,3);
t331 = pkin(2) * t201;
t336 = m(6) / 0.2e1;
t337 = m(5) / 0.2e1;
t354 = t153 * mrSges(6,3);
t53 = t373 * t301;
t54 = t89 * t354;
t359 = -t84 / 0.2e1 - t130 / 0.2e1 + t331 / 0.2e1 - t252 / 0.2e1 + t54 / 0.2e1 - t53 / 0.2e1 - t251 / 0.2e1 - t283 / 0.2e1 - (t246 + t247) * t336 - (t197 + t261) * t329 * t337 - (-t153 * t66 + t257 * t374) * mrSges(6,3) / 0.2e1 + (t373 + t374) * t301 / 0.2e1 + (-t89 - t66) * t354 / 0.2e1 - t371;
t297 = t190 * mrSges(5,3);
t209 = pkin(3) * t216 + pkin(4);
t330 = pkin(3) * t215;
t173 = t209 * t220 - t217 * t330;
t174 = t209 * t217 + t220 * t330;
t271 = pkin(3) * t337;
t226 = (-t153 * t173 + t174 * t257) * t336 + (t189 * t215 + t190 * t216) * t271;
t230 = t218 * t271 + t256 * t336;
t27 = -t226 + t230 + t358 + t259;
t272 = qJD(1) + qJD(2);
t349 = t272 * t27;
t42 = 0.2e1 * t267 + t350;
t348 = t272 * t42;
t344 = t360 * t222;
t342 = -mrSges(4,1) * t221 + mrSges(4,2) * t218;
t165 = t190 * t332;
t166 = t189 * t332;
t109 = t165 * t220 - t166 * t217;
t110 = t165 * t217 + t166 * t220;
t276 = t109 * mrSges(6,1) / 0.2e1 - t110 * mrSges(6,2) / 0.2e1;
t340 = (t109 * t173 + t110 * t174) * t336 + t165 * mrSges(5,1) / 0.2e1 - t166 * mrSges(5,2) / 0.2e1 + (t165 * t216 + t166 * t215) * t271 + t276 + (-t295 / 0.2e1 - t296 / 0.2e1) * t332;
t339 = (-t109 * t153 + t110 * t257) * mrSges(6,3) + (t165 * t190 + t166 * t189) * mrSges(5,3);
t334 = m(5) * t219;
t333 = m(6) * t219;
t318 = t27 * qJD(3) + t42 * qJD(5);
t31 = t226 + t230;
t43 = t267 - t355 / 0.2e1;
t317 = t31 * qJD(3) + t43 * qJD(5);
t316 = -t153 * t374 + t257 * t66;
t315 = -t153 * t373 + t257 * t89;
t307 = pkin(3) * qJD(3);
t298 = t189 * mrSges(5,3);
t3 = m(6) * t247 + t197 * t365 + t251 + t252 + t283 + t371;
t294 = t3 * qJD(1);
t231 = (t372 * t153 + t148) * t257 - t153 ^ 2 * Ifges(6,4);
t7 = t231 + t252;
t291 = t7 * qJD(1);
t232 = t153 * t354 + t257 * t301 + (t189 ^ 2 + t190 ^ 2) * mrSges(5,3);
t275 = t142 * t189 + t190 * t346;
t18 = m(5) * t275 + m(6) * t316 + t232;
t289 = qJD(1) * t18;
t229 = (t360 * mrSges(4,3) - mrSges(3,2)) * t222 + (-mrSges(3,1) + t154 + t91 + t342) * t219;
t15 = m(6) * (t109 * t374 + t110 * t66) + m(5) * (t142 * t166 + t165 * t346) + (t160 * t333 + t197 * t334 + m(4) * (t344 * t210 + t211 * t219) + t229) * pkin(1) + t339;
t286 = t15 * qJD(1);
t58 = -mrSges(6,1) * t174 - mrSges(6,2) * t173;
t278 = t58 * qJD(5);
t274 = t159 * t189 + t190 * t345;
t270 = t216 * t298;
t1 = t340 + t359;
t4 = t331 - t345 * t298 - t159 * t297 - t130 - t53 + t54 - t84 - t218 * t380 + (mrSges(5,3) * t159 + t309) * t190 - t261 * t365 - m(6) * t246 + (mrSges(6,3) * t373 - t148) * t257 + (-t214 + t213) * Ifges(4,4) + (t345 * mrSges(5,3) + t375 * t190 - t310) * t189 - t290 + (-mrSges(6,3) * t89 - t255 + t308) * t153;
t254 = -t1 * qJD(1) - t4 * qJD(2);
t10 = t231 + t84;
t227 = (t160 / 0.2e1 + t164 / 0.2e1) * t358 + t231;
t5 = t227 - t276;
t250 = t5 * qJD(1) + t10 * qJD(2);
t225 = (t274 + t275) * t337 + (t315 + t316) * t336 + t232;
t14 = (-m(6) / 0.2e1 - m(5) / 0.2e1) * t328 + t225;
t20 = m(5) * t274 + m(6) * t315 + t232;
t248 = -qJD(1) * t14 - qJD(2) * t20;
t233 = t58 * qJD(3);
t228 = Ifges(4,5) * t221 + Ifges(5,5) * t189 - Ifges(4,6) * t218 + Ifges(5,6) * t190 - t173 * t301 - t174 * t354 + t297 * t330 + t378;
t35 = t42 * qJD(4);
t34 = t43 * qJD(4);
t29 = t31 * qJD(4);
t26 = t27 * qJD(4);
t13 = t225 + (m(5) + m(6)) * t328 / 0.2e1;
t6 = t227 + t276;
t2 = t340 - t359;
t9 = [qJD(2) * t15 + qJD(3) * t3 + qJD(4) * t18 + qJD(5) * t7, t2 * qJD(3) + t13 * qJD(4) + t6 * qJD(5) + t286 + (0.2e1 * (t109 * t373 + t110 * t89) * t336 + 0.2e1 * (t159 * t166 + t165 * t345) * t337 + (t164 * t333 + t261 * t334 + m(4) * (-pkin(2) * t219 + t344 * pkin(7)) + t229) * pkin(1) + t339) * qJD(2), t294 + t2 * qJD(2) + (t228 + m(6) * (-t173 * t66 + t174 * t374) - t142 * mrSges(5,1) - t346 * mrSges(5,2) + t342 * t210 + t390) * qJD(3) + t29 + t388 + (-t270 + m(5) * (-t142 * t216 + t215 * t346)) * t307, qJD(2) * t13 + t289 + t317, t6 * qJD(2) + t8 * qJD(3) + t291 + t34 + t388; -qJD(3) * t1 + qJD(4) * t14 + qJD(5) * t5 - t286, -qJD(3) * t4 + qJD(4) * t20 + qJD(5) * t10, (t228 + m(6) * (-t173 * t89 + t174 * t373) - t159 * mrSges(5,1) - t345 * mrSges(5,2) + t342 * pkin(7) + t389) * qJD(3) + t29 + t387 + (-t270 + m(5) * (-t159 * t216 + t215 * t345)) * t307 + t254, -t248 + t317, t11 * qJD(3) + t250 + t34 + t387; qJD(2) * t1 - t26 - t294, -t254 - t26, t278, -t349, t233 + t278; -qJD(2) * t14 - t289 + t318, t248 + t318, t349, 0, t348; -qJD(2) * t5 - t291 - t35, -t250 - t35, -t233, -t348, 0;];
Cq = t9;
