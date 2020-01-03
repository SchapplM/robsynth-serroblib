% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:43
% EndTime: 2019-12-31 21:50:50
% DurationCPUTime: 3.30s
% Computational Cost: add. (7259->258), mult. (14488->314), div. (0->0), fcn. (13299->6), ass. (0->171)
t393 = Ifges(5,4) - Ifges(6,5);
t194 = cos(qJ(3));
t192 = sin(qJ(3));
t322 = sin(qJ(4));
t251 = t322 * t192;
t323 = cos(qJ(4));
t162 = -t194 * t323 + t251;
t252 = t323 * t192;
t163 = -t194 * t322 - t252;
t113 = mrSges(5,1) * t162 - mrSges(5,2) * t163;
t361 = -Ifges(4,4) * t192 + (Ifges(4,1) - Ifges(4,2)) * t194;
t344 = -Ifges(5,1) - Ifges(6,1);
t340 = Ifges(6,3) + t344;
t388 = t393 * t162;
t389 = t393 * t163;
t362 = -t163 * t389 + t162 * ((-Ifges(5,2) - t340) * t163 + t388);
t391 = t194 ^ 2;
t374 = t391 * Ifges(4,4);
t392 = (pkin(3) * t113 + t361) * t192 + t374 + t362;
t390 = mrSges(5,1) + mrSges(6,1);
t188 = t194 * pkin(8);
t265 = t194 * pkin(7) + t188;
t328 = -pkin(8) - pkin(7);
t338 = t328 * t251 + t265 * t323;
t355 = t338 * mrSges(6,1);
t356 = t338 * mrSges(5,1);
t124 = -t328 * t252 + t265 * t322;
t376 = t124 * mrSges(6,3);
t377 = t124 * mrSges(5,2);
t387 = -t355 - t356 + t377 - t376;
t193 = sin(qJ(2));
t318 = pkin(1) * t193;
t253 = pkin(7) + t318;
t223 = t192 * (-pkin(8) - t253);
t175 = t194 * t253;
t266 = t175 + t188;
t339 = t322 * t223 + t266 * t323;
t353 = t339 * mrSges(6,1);
t354 = t339 * mrSges(5,1);
t98 = -t323 * t223 + t266 * t322;
t378 = t98 * mrSges(6,3);
t379 = t98 * mrSges(5,2);
t386 = -t353 - t354 + t379 - t378;
t195 = cos(qJ(2));
t317 = pkin(1) * t195;
t136 = t162 * t317;
t385 = (-mrSges(6,3) / 0.2e1 + mrSges(5,2) / 0.2e1) * t136;
t384 = t355 / 0.2e1 + t356 / 0.2e1 + t376 / 0.2e1 - t377 / 0.2e1;
t383 = t353 / 0.2e1 + t354 / 0.2e1 + t378 / 0.2e1 - t379 / 0.2e1;
t135 = t163 * t317;
t171 = mrSges(4,1) * t192 + mrSges(4,2) * t194;
t261 = t322 * pkin(3);
t177 = t261 + qJ(5);
t262 = t323 * pkin(3);
t181 = -t262 - pkin(4);
t330 = m(5) * pkin(3);
t331 = m(6) / 0.2e1;
t337 = t390 * t135 / 0.2e1 + t385;
t373 = (-t135 * t181 - t136 * t177) * t331 + (t135 * t323 - t136 * t322) * t330 / 0.2e1 - t171 * t317 / 0.2e1 + t337;
t368 = -pkin(4) * t339 - qJ(5) * t98;
t367 = -t177 * t98 + t181 * t339;
t366 = t322 * t98 + t323 * t339;
t364 = -pkin(4) * t338 - qJ(5) * t124;
t109 = -pkin(4) * t163 + qJ(5) * t162;
t315 = pkin(3) * t192;
t102 = t109 + t315;
t112 = mrSges(6,1) * t162 + mrSges(6,3) * t163;
t277 = t102 * t112;
t363 = t277 + t392;
t345 = mrSges(6,2) + mrSges(5,3);
t264 = t192 ^ 2 + t391;
t185 = m(6) * qJ(5) + mrSges(6,3);
t349 = qJD(4) * t185;
t348 = t185 * qJD(5);
t304 = mrSges(4,2) * t192;
t170 = -mrSges(4,1) * t194 + t304;
t346 = -mrSges(3,1) + t112 + t113 + t170;
t236 = t261 / 0.2e1;
t285 = t162 * mrSges(6,2);
t134 = t181 * t285;
t237 = t163 * t261;
t141 = mrSges(5,3) * t237;
t238 = t162 * t262;
t336 = -t134 / 0.2e1 + t141 / 0.2e1 - mrSges(6,2) * t238 / 0.2e1;
t333 = m(5) / 0.2e1;
t332 = -m(6) / 0.2e1;
t329 = m(6) * pkin(3);
t321 = m(6) * t339;
t320 = m(6) * t338;
t319 = m(6) * t163;
t316 = pkin(2) * t171;
t284 = t163 * mrSges(6,2);
t183 = -t194 * pkin(3) - pkin(2);
t168 = t183 - t317;
t161 = t168 * t315;
t111 = -mrSges(5,1) * t163 - mrSges(5,2) * t162;
t273 = t168 * t111;
t110 = -mrSges(6,1) * t163 + mrSges(6,3) * t162;
t232 = t162 * pkin(4) + t163 * qJ(5);
t101 = t183 + t232;
t87 = t101 - t317;
t281 = t87 * t110;
t227 = t273 + t281;
t182 = -pkin(2) - t317;
t272 = t182 * t171;
t31 = t87 * t102;
t5 = m(5) * t161 + m(6) * t31 + t227 + t272 + t363;
t283 = t5 * qJD(1);
t203 = t109 * t112 + t362;
t36 = t87 * t109;
t7 = m(6) * t36 + t203 + t227;
t282 = t7 * qJD(1);
t73 = t163 * t112;
t25 = t319 * t87 + t73;
t279 = qJD(1) * t25;
t278 = t101 * t110;
t207 = t345 * (t135 * t163 + t136 * t162);
t13 = t207 + (m(6) + m(5)) * (-t135 * t98 - t136 * t339) + (m(4) * t182 + m(5) * t168 + m(6) * t87 + t346) * t318 + (-mrSges(3,2) + (m(4) * t253 + mrSges(4,3)) * t264) * t317;
t276 = t13 * qJD(1);
t271 = t183 * t111;
t257 = t284 / 0.2e1;
t267 = pkin(4) * t285 / 0.2e1 + qJ(5) * t257;
t260 = t135 * t331;
t150 = -t285 / 0.2e1;
t258 = -t284 / 0.2e1;
t235 = 0.2e1 * t150;
t234 = (Ifges(5,6) - Ifges(6,6)) * t163 + (-Ifges(6,4) - Ifges(5,5)) * t162;
t169 = t183 * t315;
t37 = t101 * t102;
t200 = t277 + (t161 + t169) * t333 + (t31 + t37) * t331;
t1 = t361 * t192 + t200 + t281 / 0.2e1 + t272 / 0.2e1 + t273 / 0.2e1 + t113 * t315 - t316 / 0.2e1 + t271 / 0.2e1 + t278 / 0.2e1 + t374 - t373 + t388 * t162 + ((-Ifges(6,3) / 0.2e1 - t340 / 0.2e1 - Ifges(5,2) - t344 / 0.2e1) * t162 - t389) * t163;
t224 = t271 + t278;
t6 = m(5) * t169 + m(6) * t37 + t224 - t316 + t363;
t231 = t1 * qJD(1) + t6 * qJD(2);
t209 = (t101 / 0.2e1 + t87 / 0.2e1) * t110 + (t183 / 0.2e1 + t168 / 0.2e1) * t111;
t40 = t101 * t109;
t196 = (t36 + t40) * t331 + t203 + t209;
t219 = m(6) * (pkin(4) * t135 - qJ(5) * t136);
t3 = -t385 - (mrSges(5,1) / 0.2e1 + mrSges(6,1) / 0.2e1) * t135 - t219 / 0.2e1 + t196;
t8 = m(6) * t40 + t203 + t224;
t230 = t3 * qJD(1) + t8 * qJD(2);
t218 = -t73 - (t101 + t87) * t319 / 0.2e1;
t18 = -t260 + t218;
t26 = t101 * t319 + t73;
t225 = qJD(1) * t18 - qJD(2) * t26;
t222 = -qJD(3) * t185 - t349;
t220 = t177 * t257 + t336 - t345 * t237 / 0.2e1;
t198 = (t366 * pkin(3) + t367) * t332 + t383;
t201 = t368 * t331 + t267 - t383;
t10 = t345 * t163 * t236 + t177 * t258 + t198 + t201 - t336;
t197 = ((t261 - t177) * t124 + (t181 + t262) * t338) * t331 + t220 - t384;
t202 = t364 * t332 + t384;
t12 = pkin(4) * t150 + qJ(5) * t258 + t197 + t202;
t206 = -t390 * t261 + (-mrSges(5,2) + mrSges(6,3)) * t262;
t63 = -(t177 * t323 + t181 * t322) * t329 - t206;
t213 = -t10 * qJD(1) + t12 * qJD(2) - t63 * qJD(3);
t167 = m(6) * t177 + mrSges(6,3);
t212 = qJD(3) * t167;
t208 = mrSges(6,2) * t232 + t234;
t204 = mrSges(5,3) * t238 + Ifges(4,5) * t194 - Ifges(4,6) * t192 + t177 * t284 - t134 + t141 + t234;
t148 = mrSges(6,3) + (0.2e1 * t236 + qJ(5)) * m(6);
t45 = t235 + t320;
t34 = t235 + t321;
t27 = t320 / 0.2e1 + t338 * t331 + t235;
t22 = t321 / 0.2e1 + t339 * t331 + t235;
t19 = -t260 - t218;
t11 = t197 - t202 + t234 + t267;
t9 = -t198 + t201 + t220 + t234;
t4 = t219 / 0.2e1 + t196 + t337;
t2 = (-pkin(2) / 0.2e1 + t182 / 0.2e1) * t171 + t200 + t209 + t373 + t392;
t14 = [qJD(2) * t13 + qJD(3) * t5 + qJD(4) * t7 + qJD(5) * t25, t2 * qJD(3) + t4 * qJD(4) + t19 * qJD(5) + t276 + (t207 + 0.2e1 * (t333 + t331) * (-t124 * t135 - t136 * t338) + ((-m(4) * pkin(2) + m(5) * t183 + m(6) * t101 + t346) * t193 + (-mrSges(3,2) + (m(4) * pkin(7) + mrSges(4,3)) * t264) * t195) * pkin(1)) * qJD(2), t283 + t2 * qJD(2) + (m(6) * t367 - mrSges(4,1) * t175 + t253 * t304 - t366 * t330 + t204 + t386) * qJD(3) + t9 * qJD(4) + t22 * qJD(5), t282 + t4 * qJD(2) + t9 * qJD(3) + (m(6) * t368 + t208 + t386) * qJD(4) + t34 * qJD(5), qJD(2) * t19 + qJD(3) * t22 + qJD(4) * t34 + t279; qJD(3) * t1 + qJD(4) * t3 - qJD(5) * t18 - t276, qJD(3) * t6 + qJD(4) * t8 + qJD(5) * t26, (m(6) * (-t124 * t177 + t181 * t338) + t204 + (-t124 * t322 - t323 * t338) * t330 + t170 * pkin(7) + t387) * qJD(3) + t11 * qJD(4) + t27 * qJD(5) + t231, t11 * qJD(3) + (m(6) * t364 + t208 + t387) * qJD(4) + t45 * qJD(5) + t230, qJD(3) * t27 + qJD(4) * t45 - t225; -qJD(2) * t1 - qJD(4) * t10 - t283, qJD(4) * t12 - t231, -qJD(4) * t63 + qJD(5) * t167, ((-pkin(4) * t322 + qJ(5) * t323) * t329 + t206) * qJD(4) + t148 * qJD(5) + t213, qJD(4) * t148 + t212; -qJD(2) * t3 + qJD(3) * t10 - t282, -qJD(3) * t12 - t230, -t213 + t348, t348, -t222; qJD(2) * t18 - t279, t225, -t212 - t349, t222, 0;];
Cq = t14;
