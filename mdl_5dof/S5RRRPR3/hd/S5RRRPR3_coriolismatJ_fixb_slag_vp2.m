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
% m [6x1]
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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:42:31
% EndTime: 2022-01-20 11:42:39
% DurationCPUTime: 3.81s
% Computational Cost: add. (11144->250), mult. (21423->331), div. (0->0), fcn. (22873->8), ass. (0->160)
t228 = sin(qJ(3));
t229 = sin(qJ(2));
t329 = t229 * pkin(1);
t218 = pkin(7) + t329;
t272 = qJ(4) + t218;
t186 = t272 * t228;
t231 = cos(qJ(3));
t187 = t272 * t231;
t225 = sin(pkin(9));
t226 = cos(pkin(9));
t139 = -t225 * t186 + t226 * t187;
t193 = -t225 * t228 + t226 * t231;
t184 = t193 * pkin(8);
t104 = t184 + t139;
t227 = sin(qJ(5));
t230 = cos(qJ(5));
t194 = -t225 * t231 - t226 * t228;
t185 = t194 * pkin(8);
t351 = -t226 * t186 - t225 * t187;
t363 = t351 + t185;
t390 = -t104 * t227 + t230 * t363;
t62 = t104 * t230 + t227 * t363;
t404 = -t62 * mrSges(6,1) - t390 * mrSges(6,2);
t320 = -qJ(4) - pkin(7);
t202 = t320 * t228;
t204 = t320 * t231;
t161 = t225 * t202 - t226 * t204;
t123 = t184 + t161;
t350 = t226 * t202 + t225 * t204;
t364 = t350 + t185;
t389 = -t123 * t227 + t230 * t364;
t84 = t123 * t230 + t227 * t364;
t403 = -t84 * mrSges(6,1) - t389 * mrSges(6,2);
t150 = t193 * t227 - t194 * t230;
t251 = t230 * t193 + t194 * t227;
t392 = Ifges(6,5) * t251 - Ifges(6,6) * t150;
t8 = t392 + t404;
t402 = t8 * qJD(5);
t11 = t392 + t403;
t401 = t11 * qJD(5);
t330 = pkin(3) * t228;
t169 = -pkin(4) * t194 + t330;
t400 = m(6) * t169;
t338 = t251 / 0.2e1;
t395 = 0.2e1 * t338;
t365 = t150 * mrSges(6,3);
t388 = t193 / 0.2e1;
t387 = t231 / 0.2e1;
t220 = -t231 * pkin(3) - pkin(2);
t166 = -t193 * pkin(4) + t220;
t232 = cos(qJ(2));
t333 = pkin(1) * t232;
t162 = t166 - t333;
t366 = t150 * mrSges(6,1);
t372 = t251 * mrSges(6,2) + t366;
t297 = t162 * t372;
t296 = t166 * t372;
t383 = (t162 / 0.2e1 + t166 / 0.2e1) * t372;
t381 = Ifges(6,4) * t251 * t395 + (Ifges(6,1) * t395 - Ifges(6,4) * t150 + (-t251 / 0.2e1 - t338) * Ifges(6,2)) * t150;
t382 = (Ifges(5,2) * t388 - Ifges(5,4) * t194 - Ifges(5,1) * t193 / 0.2e1) * t194 + t381;
t260 = t366 / 0.2e1;
t378 = m(5) * t330;
t373 = t228 ^ 2 + t231 ^ 2;
t368 = -Ifges(4,2) * t231 / 0.2e1 - Ifges(4,4) * t228 + Ifges(4,1) * t387;
t151 = -t194 * mrSges(5,1) + t193 * mrSges(5,2);
t244 = m(5) * (t193 * t225 + t194 * t226) * pkin(3);
t343 = m(5) / 0.2e1;
t252 = t330 * t343;
t217 = pkin(3) * t226 + pkin(4);
t331 = pkin(3) * t225;
t177 = t217 * t227 + t230 * t331;
t281 = t177 * t251;
t176 = t217 * t230 - t227 * t331;
t282 = t176 * t150;
t26 = t252 + 0.2e1 * (t169 / 0.4e1 + t282 / 0.4e1 - t281 / 0.4e1) * m(6) - t244 / 0.2e1 + t372 + t151;
t267 = qJD(1) + qJD(2);
t355 = t267 * t26;
t41 = mrSges(6,2) * t395 + 0.2e1 * t260;
t354 = t267 * t41;
t349 = t373 * t232;
t348 = -mrSges(4,1) * t231 + mrSges(4,2) * t228;
t167 = t194 * t333;
t168 = t193 * t333;
t109 = t167 * t230 - t168 * t227;
t110 = t167 * t227 + t168 * t230;
t271 = t109 * mrSges(6,1) / 0.2e1 - t110 * mrSges(6,2) / 0.2e1;
t346 = (-t109 * t150 + t110 * t251) * mrSges(6,3) + (t167 * t194 + t168 * t193) * mrSges(5,3);
t183 = Ifges(5,4) * t193;
t153 = Ifges(5,2) * t194 + t183;
t156 = -Ifges(5,1) * t194 + t183;
t222 = Ifges(4,4) * t231;
t207 = -Ifges(4,2) * t228 + t222;
t208 = Ifges(4,1) * t228 + t222;
t345 = (t153 + t156) * t388 + (t207 + t208) * t387 + t368 * t228 + t382;
t342 = -m(6) / 0.2e1;
t341 = m(6) / 0.2e1;
t335 = m(5) * t229;
t334 = m(6) * t229;
t292 = t231 * mrSges(4,2);
t293 = t228 * mrSges(4,1);
t205 = t292 + t293;
t332 = pkin(2) * t205;
t319 = t26 * qJD(3) + t41 * qJD(5);
t30 = t244 / 0.2e1 + t252 + (t281 - t282 + t169) * t341;
t42 = t260 - t366 / 0.2e1;
t318 = t30 * qJD(3) + t42 * qJD(5);
t317 = -t150 * t390 + t251 * t62;
t316 = -t150 * t389 + t251 * t84;
t310 = pkin(3) * qJD(3);
t302 = t251 * mrSges(6,3);
t201 = t220 - t333;
t152 = -mrSges(5,1) * t193 - mrSges(5,2) * t194;
t87 = -mrSges(6,1) * t251 + mrSges(6,2) * t150;
t289 = t152 * t330 + t169 * t87;
t238 = t289 + t345;
t219 = -pkin(2) - t333;
t279 = t219 * t205;
t280 = t201 * t151;
t3 = t162 * t400 + t201 * t378 + t238 + t279 + t280 + t297;
t291 = t3 * qJD(1);
t7 = t381 + t297;
t290 = t7 * qJD(1);
t242 = t150 * t365 + t251 * t302 + (t193 ^ 2 + t194 ^ 2) * mrSges(5,3);
t270 = t139 * t193 + t194 * t351;
t18 = m(5) * t270 + m(6) * t317 + t242;
t288 = qJD(1) * t18;
t240 = (t373 * mrSges(4,3) - mrSges(3,2)) * t232 + (-mrSges(3,1) + t152 + t87 + t348) * t229;
t15 = m(6) * (t109 * t390 + t110 * t62) + m(5) * (t139 * t168 + t167 * t351) + (t162 * t334 + t201 * t335 + m(4) * (t349 * t218 + t219 * t229) + t240) * pkin(1) + t346;
t285 = t15 * qJD(1);
t278 = t220 * t151;
t54 = -mrSges(6,1) * t177 - mrSges(6,2) * t176;
t273 = t54 * qJD(5);
t269 = t161 * t193 + t194 * t350;
t266 = t226 * t193 * mrSges(5,3);
t265 = -t333 / 0.2e1;
t253 = (t201 + t220) * t228;
t236 = (t162 + t166) * t169 * t342 - t289;
t237 = (t109 * t176 + t110 * t177) * t341 + t167 * mrSges(5,1) / 0.2e1 - t168 * mrSges(5,2) / 0.2e1 + t271;
t245 = m(5) * (t167 * t226 + t168 * t225);
t2 = (-t201 / 0.2e1 - t220 / 0.2e1) * t151 + (mrSges(4,1) * t265 - t368) * t228 + (mrSges(4,2) * t265 - t208 / 0.2e1 - t207 / 0.2e1) * t231 + 0.2e1 * (t245 / 0.4e1 - m(5) * t253 / 0.4e1) * pkin(3) + t237 + t236 + (pkin(2) / 0.2e1 - t219 / 0.2e1) * t205 + (-t156 / 0.2e1 - t153 / 0.2e1) * t193 - t383 - t382;
t4 = t166 * t400 + t220 * t378 + t238 + t278 + t296 - t332;
t250 = -t2 * qJD(1) + t4 * qJD(2);
t10 = t381 + t296;
t233 = t381 + t383;
t5 = t233 - t271;
t249 = t5 * qJD(1) + t10 * qJD(2);
t234 = (t269 + t270) * t343 + (t316 + t317) * t341 + t242;
t14 = (t342 - m(5) / 0.2e1) * t329 + t234;
t20 = m(5) * t269 + m(6) * t316 + t242;
t247 = -qJD(1) * t14 - qJD(2) * t20;
t243 = t54 * qJD(3);
t239 = Ifges(4,5) * t231 + Ifges(5,5) * t193 - Ifges(4,6) * t228 - t176 * t302 - t177 * t365 + t392 + (mrSges(5,3) * t331 + Ifges(5,6)) * t194;
t34 = t41 * qJD(4);
t33 = t42 * qJD(4);
t28 = t30 * qJD(4);
t25 = t26 * qJD(4);
t13 = t234 + (m(5) + m(6)) * t329 / 0.2e1;
t6 = t233 + t271;
t1 = t296 / 0.2e1 + t297 / 0.2e1 + (-t292 / 0.2e1 - t293 / 0.2e1) * t333 + t237 + t345 - t236 + t280 / 0.2e1 + t278 / 0.2e1 + t279 / 0.2e1 - t332 / 0.2e1 + (t253 * t343 + t245 / 0.2e1) * pkin(3) + (mrSges(6,3) * t338 - t302 / 0.2e1) * (t389 + t390);
t9 = [qJD(2) * t15 + qJD(3) * t3 + qJD(4) * t18 + qJD(5) * t7, t1 * qJD(3) + t13 * qJD(4) + t6 * qJD(5) + t285 + (0.2e1 * (t109 * t389 + t110 * t84) * t341 + 0.2e1 * (t161 * t168 + t167 * t350) * t343 + (t166 * t334 + t220 * t335 + m(4) * (-pkin(2) * t229 + t349 * pkin(7)) + t240) * pkin(1) + t346) * qJD(2), t291 + t1 * qJD(2) + (t239 + m(6) * (-t176 * t62 + t177 * t390) - t139 * mrSges(5,1) - t351 * mrSges(5,2) + t348 * t218 + t404) * qJD(3) + t28 + t402 + (-t266 + m(5) * (-t139 * t226 + t225 * t351)) * t310, qJD(2) * t13 + t288 + t318, t6 * qJD(2) + t8 * qJD(3) + t290 + t33 + t402; -qJD(3) * t2 + qJD(4) * t14 + qJD(5) * t5 - t285, qJD(3) * t4 + qJD(4) * t20 + qJD(5) * t10, (t239 + m(6) * (-t176 * t84 + t177 * t389) - t161 * mrSges(5,1) - t350 * mrSges(5,2) + t348 * pkin(7) + t403) * qJD(3) + t28 + t401 + (-t266 + m(5) * (-t161 * t226 + t225 * t350)) * t310 + t250, -t247 + t318, t11 * qJD(3) + t249 + t33 + t401; qJD(2) * t2 - t25 - t291, -t25 - t250, t273, -t355, t243 + t273; -qJD(2) * t14 - t288 + t319, t247 + t319, t355, 0, t354; -qJD(2) * t5 - t290 - t34, -t249 - t34, -t243, -t354, 0;];
Cq = t9;
