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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:42:23
% EndTime: 2019-12-05 18:42:31
% DurationCPUTime: 3.89s
% Computational Cost: add. (11144->234), mult. (21423->300), div. (0->0), fcn. (22873->8), ass. (0->151)
t219 = sin(qJ(3));
t220 = sin(qJ(2));
t316 = t220 * pkin(1);
t210 = pkin(7) + t316;
t267 = qJ(4) + t210;
t180 = t267 * t219;
t222 = cos(qJ(3));
t181 = t267 * t222;
t217 = sin(pkin(9));
t281 = cos(pkin(9));
t138 = -t217 * t180 + t281 * t181;
t188 = -t217 * t219 + t281 * t222;
t178 = t188 * pkin(8);
t101 = t178 + t138;
t218 = sin(qJ(5));
t221 = cos(qJ(5));
t189 = -t217 * t222 - t281 * t219;
t179 = t189 * pkin(8);
t345 = -t281 * t180 - t217 * t181;
t357 = t345 + t179;
t391 = -t101 * t218 + t221 * t357;
t62 = t101 * t221 + t218 * t357;
t404 = -t62 * mrSges(6,1) - t391 * mrSges(6,2);
t307 = -qJ(4) - pkin(7);
t198 = t307 * t219;
t200 = t307 * t222;
t155 = t217 * t198 - t281 * t200;
t121 = t178 + t155;
t344 = t281 * t198 + t217 * t200;
t358 = t344 + t179;
t390 = -t121 * t218 + t221 * t358;
t83 = t121 * t221 + t218 * t358;
t403 = -t83 * mrSges(6,1) - t390 * mrSges(6,2);
t148 = t188 * t218 - t189 * t221;
t245 = t221 * t188 + t189 * t218;
t394 = Ifges(6,5) * t245 - Ifges(6,6) * t148;
t8 = t394 + t404;
t402 = t8 * qJD(5);
t11 = t394 + t403;
t401 = t11 * qJD(5);
t363 = t148 * mrSges(6,3);
t353 = t245 / 0.2e1;
t392 = 0.2e1 * t353;
t389 = t222 ^ 2;
t212 = -t222 * pkin(3) - pkin(2);
t160 = -t188 * pkin(4) + t212;
t223 = cos(qJ(2));
t320 = pkin(1) * t223;
t156 = t160 - t320;
t364 = t148 * mrSges(6,1);
t374 = t245 * mrSges(6,2) + t364;
t287 = t156 * t374;
t286 = t160 * t374;
t386 = (t156 / 0.2e1 + t160 / 0.2e1) * t374;
t150 = -mrSges(5,1) * t188 - mrSges(5,2) * t189;
t177 = Ifges(5,4) * t188;
t297 = Ifges(5,4) * t189;
t366 = -Ifges(5,1) + Ifges(5,2);
t370 = -Ifges(4,4) * t219 + (Ifges(4,1) - Ifges(4,2)) * t222;
t379 = t389 * Ifges(4,4);
t385 = (pkin(3) * t150 + t370) * t219 + t379 + (t366 * t189 + t177) * t188 - t297 * t189;
t255 = t364 / 0.2e1;
t377 = t219 ^ 2 + t389;
t161 = t189 * t320;
t162 = t188 * t320;
t106 = t161 * t221 - t162 * t218;
t107 = t161 * t218 + t162 * t221;
t250 = t281 * pkin(3);
t209 = t250 + pkin(4);
t318 = pkin(3) * t217;
t170 = t209 * t221 - t218 * t318;
t171 = t209 * t218 + t221 * t318;
t201 = mrSges(4,1) * t219 + mrSges(4,2) * t222;
t334 = m(5) * pkin(3);
t261 = t334 / 0.2e1;
t266 = t106 * mrSges(6,1) / 0.2e1 - t107 * mrSges(6,2) / 0.2e1;
t335 = m(6) / 0.2e1;
t376 = (t106 * t170 + t107 * t171) * t335 + t161 * mrSges(5,1) / 0.2e1 - t162 * mrSges(5,2) / 0.2e1 + (t281 * t161 + t162 * t217) * t261 - t201 * t320 / 0.2e1 + t266;
t317 = pkin(3) * t219;
t163 = -pkin(4) * t189 + t317;
t230 = Ifges(6,4) * t245 * t392 + (-Ifges(6,4) * t148 + (Ifges(6,1) - Ifges(6,2)) * t353 + (-Ifges(6,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t245) * t148;
t86 = -mrSges(6,1) * t245 + mrSges(6,2) * t148;
t228 = t163 * t86 + t230;
t371 = t228 + t385;
t149 = -t189 * mrSges(5,1) + t188 * mrSges(5,2);
t232 = (-t148 * t170 + t171 * t245) * t335 + (t188 * t217 + t281 * t189) * t261;
t237 = t163 * t335 + t219 * t261;
t26 = t149 - t232 + t237 + t374;
t262 = qJD(1) + qJD(2);
t349 = t262 * t26;
t41 = mrSges(6,2) * t392 + 0.2e1 * t255;
t348 = t262 * t41;
t342 = t377 * t223;
t341 = -mrSges(4,1) * t222 + mrSges(4,2) * t219;
t337 = (t161 * t189 + t162 * t188) * mrSges(5,3) + (-t106 * t148 + t107 * t245) * mrSges(6,3);
t336 = m(5) / 0.2e1;
t322 = m(5) * t220;
t321 = m(6) * t220;
t319 = pkin(2) * t201;
t306 = t26 * qJD(3) + t41 * qJD(5);
t30 = t232 + t237;
t42 = t255 - t364 / 0.2e1;
t305 = t30 * qJD(3) + t42 * qJD(5);
t304 = -t148 * t391 + t245 * t62;
t303 = -t148 * t390 + t245 * t83;
t290 = t245 * mrSges(6,3);
t112 = t156 * t163;
t196 = t212 - t320;
t186 = t196 * t317;
t211 = -pkin(2) - t320;
t273 = t211 * t201;
t274 = t196 * t149;
t3 = m(5) * t186 + m(6) * t112 + t273 + t274 + t287 + t371;
t283 = t3 * qJD(1);
t7 = t230 + t287;
t282 = t7 * qJD(1);
t235 = t148 * t363 + t245 * t290 + (t188 ^ 2 + t189 ^ 2) * mrSges(5,3);
t265 = t138 * t188 + t189 * t345;
t18 = m(5) * t265 + m(6) * t304 + t235;
t280 = qJD(1) * t18;
t233 = (t377 * mrSges(4,3) - mrSges(3,2)) * t223 + (-mrSges(3,1) + t150 + t86 + t341) * t220;
t15 = m(6) * (t106 * t391 + t107 * t62) + m(5) * (t138 * t162 + t161 * t345) + (t156 * t321 + t196 * t322 + m(4) * (t342 * t210 + t211 * t220) + t233) * pkin(1) + t337;
t277 = t15 * qJD(1);
t272 = t212 * t149;
t54 = -mrSges(6,1) * t171 - mrSges(6,2) * t170;
t268 = t54 * qJD(5);
t264 = t155 * t188 + t189 * t344;
t123 = t160 * t163;
t197 = t212 * t317;
t224 = (t186 + t197) * t336 + (t112 + t123) * t335 + t228;
t2 = t224 + t150 * t317 + t274 / 0.2e1 + t272 / 0.2e1 + t273 / 0.2e1 + t286 / 0.2e1 + t287 / 0.2e1 + t188 * t177 + t379 - t319 / 0.2e1 + t370 * t219 + (t188 * t366 - t297) * t189 - t376 + (mrSges(6,3) * t353 - t290 / 0.2e1) * (t390 + t391);
t4 = m(5) * t197 + m(6) * t123 + t272 + t286 - t319 + t371;
t243 = t2 * qJD(1) + t4 * qJD(2);
t10 = t230 + t286;
t227 = t230 + t386;
t6 = t227 - t266;
t242 = t6 * qJD(1) + t10 * qJD(2);
t229 = (t264 + t265) * t336 + (t303 + t304) * t335 + t235;
t14 = (-m(6) / 0.2e1 - m(5) / 0.2e1) * t316 + t229;
t20 = m(5) * t264 + m(6) * t303 + t235;
t240 = -qJD(1) * t14 - qJD(2) * t20;
t236 = t54 * qJD(3);
t231 = Ifges(4,5) * t222 - Ifges(4,6) * t219 - t170 * t290 - t171 * t363 + t394 + (-mrSges(5,3) * t250 + Ifges(5,5)) * t188 + (mrSges(5,3) * t318 + Ifges(5,6)) * t189;
t34 = t41 * qJD(4);
t33 = t42 * qJD(4);
t28 = t30 * qJD(4);
t25 = t26 * qJD(4);
t13 = t229 + (m(5) + m(6)) * t316 / 0.2e1;
t5 = t227 + t266;
t1 = t224 + t386 + (t212 / 0.2e1 + t196 / 0.2e1) * t149 + (-pkin(2) / 0.2e1 + t211 / 0.2e1) * t201 + t376 + t385;
t9 = [qJD(2) * t15 + qJD(3) * t3 + qJD(4) * t18 + qJD(5) * t7, t1 * qJD(3) + t13 * qJD(4) + t5 * qJD(5) + t277 + (0.2e1 * (t106 * t390 + t107 * t83) * t335 + 0.2e1 * (t155 * t162 + t161 * t344) * t336 + (t160 * t321 + t212 * t322 + m(4) * (-pkin(2) * t220 + t342 * pkin(7)) + t233) * pkin(1) + t337) * qJD(2), t283 + t1 * qJD(2) + (t231 + m(6) * (-t170 * t62 + t171 * t391) + (-t138 * t281 + t217 * t345) * t334 - t138 * mrSges(5,1) - t345 * mrSges(5,2) + t341 * t210 + t404) * qJD(3) + t28 + t402, qJD(2) * t13 + t280 + t305, t5 * qJD(2) + t8 * qJD(3) + t282 + t33 + t402; qJD(3) * t2 + qJD(4) * t14 + qJD(5) * t6 - t277, qJD(3) * t4 + qJD(4) * t20 + qJD(5) * t10, (t231 + m(6) * (-t170 * t83 + t171 * t390) + (-t155 * t281 + t217 * t344) * t334 - t344 * mrSges(5,2) - t155 * mrSges(5,1) + t341 * pkin(7) + t403) * qJD(3) + t28 + t401 + t243, -t240 + t305, t11 * qJD(3) + t242 + t33 + t401; -qJD(2) * t2 - t25 - t283, -t243 - t25, t268, -t349, t236 + t268; -qJD(2) * t14 - t280 + t306, t240 + t306, t349, 0, t348; -qJD(2) * t6 - t282 - t34, -t242 - t34, -t236, -t348, 0;];
Cq = t9;
