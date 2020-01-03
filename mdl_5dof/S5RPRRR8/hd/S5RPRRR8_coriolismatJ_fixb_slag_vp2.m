% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:05:45
% EndTime: 2019-12-31 19:05:50
% DurationCPUTime: 2.47s
% Computational Cost: add. (6622->222), mult. (12121->287), div. (0->0), fcn. (11242->6), ass. (0->145)
t205 = sin(qJ(5));
t206 = sin(qJ(4));
t208 = cos(qJ(4));
t325 = cos(qJ(5));
t177 = t205 * t208 + t325 * t206;
t209 = cos(qJ(3));
t148 = t177 * t209;
t229 = -t205 * t206 + t325 * t208;
t149 = t229 * t209;
t334 = m(6) * pkin(4);
t263 = t334 / 0.2e1;
t268 = t148 * mrSges(6,1) / 0.2e1 + t149 * mrSges(6,2) / 0.2e1;
t289 = t208 * mrSges(5,2);
t290 = t206 * mrSges(5,1);
t340 = -t290 / 0.2e1 - t289 / 0.2e1;
t216 = (-t325 * t148 + t149 * t205) * t263 - t268 + t340 * t209;
t186 = t289 + t290;
t272 = t209 * t186;
t322 = pkin(4) * t206;
t357 = m(6) * t322;
t358 = t209 * t357;
t111 = mrSges(6,1) * t177 + mrSges(6,2) * t229;
t350 = t209 * t111;
t373 = -t350 / 0.2e1;
t361 = t373 - t358 / 0.2e1;
t375 = -t272 / 0.2e1 + t216 + t361;
t376 = t375 * qJD(4);
t112 = -mrSges(6,1) * t229 + mrSges(6,2) * t177;
t359 = -t208 / 0.2e1;
t339 = Ifges(5,2) * t359 - Ifges(5,4) * t206 + Ifges(5,1) * t208 / 0.2e1;
t374 = t112 * t322 + t339 * t206;
t305 = Ifges(6,4) * t229;
t115 = -Ifges(6,2) * t177 + t305;
t304 = Ifges(6,4) * t177;
t116 = Ifges(6,2) * t229 + t304;
t118 = -Ifges(6,1) * t229 + t304;
t120 = Ifges(6,1) * t177 + t305;
t360 = -(t116 / 0.2e1 + t118 / 0.2e1) * t177 + (t120 / 0.2e1 + t115 / 0.2e1) * t229;
t207 = sin(qJ(3));
t210 = -pkin(1) - pkin(2);
t184 = t209 * qJ(2) + t207 * t210;
t181 = -pkin(7) + t184;
t314 = pkin(8) - t181;
t140 = t314 * t206;
t141 = t314 * t208;
t238 = t325 * t140 + t141 * t205;
t266 = -Ifges(6,5) * t229 + Ifges(6,6) * t177;
t63 = t205 * t140 - t325 * t141;
t19 = -t63 * mrSges(6,1) - t238 * mrSges(6,2) + t266;
t371 = t19 * qJD(5);
t333 = -pkin(8) - pkin(7);
t191 = t333 * t206;
t192 = t333 * t208;
t129 = t205 * t191 - t325 * t192;
t237 = t325 * t191 + t192 * t205;
t28 = -t129 * mrSges(6,1) - t237 * mrSges(6,2) - t266;
t370 = t28 * qJD(5);
t306 = Ifges(5,4) * t208;
t188 = -Ifges(5,2) * t206 + t306;
t189 = Ifges(5,1) * t206 + t306;
t343 = (t189 / 0.2e1 + t188 / 0.2e1) * t208;
t367 = t343 + t374;
t281 = t149 * t229;
t109 = mrSges(6,3) * t281;
t185 = -mrSges(5,1) * t208 + mrSges(5,2) * t206;
t174 = t207 * t185;
t351 = t207 * t112;
t366 = t174 + t109 + t351;
t146 = t229 * t207;
t147 = t177 * t207;
t241 = -t146 * mrSges(6,1) + t147 * mrSges(6,2);
t365 = t241 * qJD(5);
t362 = -t174 / 0.2e1 - t351 / 0.2e1;
t183 = -t207 * qJ(2) + t209 * t210;
t180 = pkin(3) - t183;
t320 = t208 * pkin(4);
t152 = t180 + t320;
t353 = t152 * t111;
t197 = -pkin(3) - t320;
t352 = t197 * t111;
t335 = m(6) / 0.2e1;
t349 = (t152 - t197) * t335;
t225 = t183 * t229;
t292 = t177 * mrSges(6,3);
t293 = t229 * mrSges(6,3);
t264 = t206 ^ 2 + t208 ^ 2;
t337 = t264 * mrSges(5,3);
t93 = t177 * t183;
t346 = t225 * t293 + t93 * t292 + (-mrSges(4,2) + t337) * t183 + (-mrSges(4,1) + t112 + t185) * t184;
t338 = -mrSges(6,2) * t225 / 0.2e1 - t93 * mrSges(6,1) / 0.2e1;
t336 = m(5) / 0.2e1;
t275 = t207 * t209;
t332 = m(6) * (t146 * t149 + t147 * t148 - t275);
t324 = m(5) * (-0.1e1 + t264) * t275;
t323 = pkin(3) * t186;
t321 = t205 * pkin(4);
t10 = -t353 + t360;
t279 = t180 * t186;
t5 = -t152 * t357 + t10 - t279 + t367;
t287 = t5 * qJD(1);
t240 = t264 * t183;
t9 = -m(6) * (t152 * t184 + t63 * t225 - t238 * t93) - m(5) * (t180 * t184 + t181 * t240) + t346;
t286 = t9 * qJD(1);
t285 = t10 * qJD(1);
t284 = t375 * qJD(1);
t253 = t350 / 0.2e1;
t13 = t253 + t268;
t283 = t13 * qJD(1);
t282 = t148 * t177;
t227 = -t148 * t238 + t63 * t149 + t152 * t207;
t233 = -t207 * mrSges(4,1) - t209 * mrSges(4,2) + mrSges(6,3) * t282;
t239 = t264 * t209;
t278 = t184 * t209;
t18 = -mrSges(3,3) - m(3) * qJ(2) + mrSges(5,3) * t239 - m(6) * t227 - m(5) * (t180 * t207 + t181 * t239) - m(4) * (-t183 * t207 + t278) + t233 + t366;
t280 = t18 * qJD(1);
t257 = pkin(4) * t325;
t245 = t281 / 0.2e1;
t34 = t324 + t332;
t211 = -t109 / 0.2e1 + ((t264 * t181 - t184) * t209 + (t180 + t240) * t207) * t336 + (t146 * t225 + t147 * t93 + t227 - t278) * t335 + t362;
t219 = t209 * t337 + t233;
t232 = (-pkin(3) * t207 + pkin(7) * t239) * t336 + (t129 * t149 - t148 * t237 + t197 * t207) * t335;
t6 = mrSges(6,3) * t245 - t211 + t219 + t232 - t362;
t234 = -t6 * qJD(1) + t34 * qJD(2);
t231 = t353 / 0.2e1 - t352 / 0.2e1;
t213 = (t205 * t225 - t325 * t93) * t263 + t340 * t183;
t3 = t231 - (t115 + t120) * t229 / 0.2e1 + (t116 + t118) * t177 / 0.2e1 - t338;
t1 = t3 - t213 + t279 / 0.2e1 + t323 / 0.2e1 + t322 * t349 + (t189 + t188) * t359 - t374;
t214 = t272 / 0.2e1 + t216;
t15 = t214 - t361;
t17 = t352 + t360;
t8 = t197 * t357 + t17 - t323 + t367;
t228 = t1 * qJD(1) - t15 * qJD(2) + t8 * qJD(3);
t22 = t373 + t268;
t226 = t3 * qJD(1) + t22 * qJD(2) + t17 * qJD(3);
t182 = (mrSges(6,1) * t205 + mrSges(6,2) * t325) * pkin(4);
t221 = t182 * qJD(4);
t217 = Ifges(5,5) * t208 - Ifges(5,6) * t206 - t257 * t293 - t292 * t321;
t4 = t231 + t338 - t360;
t173 = t182 * qJD(5);
t23 = t373 - t268;
t14 = t253 - t268;
t12 = t358 / 0.2e1 + t214 + t253;
t7 = -t148 * t292 / 0.2e1 + (t112 / 0.2e1 + t185 / 0.2e1) * t207 + (t245 + t282 / 0.2e1) * mrSges(6,3) + t211 + t232;
t2 = -t343 + (t180 / 0.2e1 + pkin(3) / 0.2e1) * t186 + ((-t112 + t349) * pkin(4) - t339) * t206 + t4 + t213;
t11 = [-qJD(2) * t18 - qJD(3) * t9 + qJD(4) * t5 + qJD(5) * t10, -t280 + 0.2e1 * (t332 / 0.2e1 + t324 / 0.2e1) * qJD(2) + t7 * qJD(3) + t12 * qJD(4) + t14 * qJD(5), -t286 + t7 * qJD(2) + (m(6) * (t129 * t225 + t197 * t184 - t237 * t93) + m(5) * (-pkin(3) * t184 + pkin(7) * t240) + t346) * qJD(3) + t2 * qJD(4) + t4 * qJD(5), t287 + t12 * qJD(2) + t2 * qJD(3) + ((t205 * t238 - t325 * t63) * t334 - t217 + t185 * t181 + t19) * qJD(4) + t371, t14 * qJD(2) + t4 * qJD(3) + t19 * qJD(4) + t285 + t371; -qJD(3) * t6 + qJD(5) * t13 + t280 - t376, t34 * qJD(3), t376 + t23 * qJD(5) + t234 + (t219 + 0.2e1 * t232 + t366) * qJD(3), -t284 + t375 * qJD(3) + (m(6) * (-t146 * t257 - t147 * t321) + t174 + t241) * qJD(4) + t365, t23 * qJD(3) + qJD(4) * t241 + t283 + t365; qJD(2) * t6 + qJD(4) * t1 + qJD(5) * t3 + t286, -qJD(4) * t15 + qJD(5) * t22 - t234, qJD(4) * t8 + qJD(5) * t17, ((-t129 * t325 + t205 * t237) * t334 + t217 + t185 * pkin(7) + t28) * qJD(4) + t370 + t228, t28 * qJD(4) + t226 + t370; qJD(2) * t375 - qJD(3) * t1 - t287, t15 * qJD(3) + t284, -t228, -t173, -t173 - t221; -qJD(2) * t13 - qJD(3) * t3 - t285, -t22 * qJD(3) - t283, -t226, t221, 0;];
Cq = t11;
