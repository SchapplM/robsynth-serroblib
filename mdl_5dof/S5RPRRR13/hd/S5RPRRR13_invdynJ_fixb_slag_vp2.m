% Calculate vector of inverse dynamics joint torques for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR13_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR13_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:33
% EndTime: 2019-12-31 19:14:57
% DurationCPUTime: 12.92s
% Computational Cost: add. (4701->561), mult. (9362->782), div. (0->0), fcn. (5722->10), ass. (0->261)
t183 = cos(qJ(4));
t166 = pkin(4) * t183 + pkin(3);
t361 = m(6) * t166;
t186 = -pkin(8) - pkin(7);
t360 = -m(6) * t186 + mrSges(5,3) + mrSges(6,3);
t184 = cos(qJ(3));
t359 = t184 / 0.2e1;
t179 = sin(qJ(4));
t180 = sin(qJ(3));
t170 = t180 * qJD(1);
t238 = t179 * t170;
t239 = qJD(4) * t186;
t224 = pkin(3) * t184 + pkin(7) * t180;
t145 = t224 * qJD(1);
t187 = -pkin(1) - pkin(6);
t162 = qJD(1) * t187 + qJD(2);
t263 = t183 * t184;
t76 = t145 * t179 + t162 * t263;
t358 = -pkin(8) * t238 + t179 * t239 - t76;
t267 = t180 * t183;
t247 = pkin(8) * t267;
t270 = t179 * t184;
t75 = t145 * t183 - t162 * t270;
t357 = t183 * t239 - (pkin(4) * t184 + t247) * qJD(1) - t75;
t254 = qJD(4) * t179;
t356 = t254 + t238;
t181 = sin(qJ(1));
t185 = cos(qJ(1));
t355 = g(1) * t181 - g(2) * t185;
t324 = m(6) * pkin(4);
t354 = m(5) + m(6);
t249 = qJD(1) * qJD(3);
t149 = qJDD(1) * t184 - t180 * t249;
t257 = qJD(3) * t183;
t259 = qJD(1) * t184;
t138 = -t179 * t259 + t257;
t63 = qJD(4) * t138 + qJDD(3) * t179 + t149 * t183;
t139 = qJD(3) * t179 + t183 * t259;
t64 = -qJD(4) * t139 + qJDD(3) * t183 - t149 * t179;
t30 = -mrSges(5,1) * t64 + mrSges(5,2) * t63;
t353 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t149 - t30;
t352 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t138 + mrSges(5,2) * t139 + mrSges(4,3) * t259;
t178 = sin(qJ(5));
t182 = cos(qJ(5));
t205 = t178 * t179 - t182 * t183;
t330 = qJD(4) + qJD(5);
t189 = t330 * t205;
t287 = Ifges(4,4) * t180;
t217 = t184 * Ifges(4,1) - t287;
t131 = Ifges(5,4) * t138;
t164 = t170 + qJD(4);
t55 = Ifges(5,1) * t139 + Ifges(5,5) * t164 + t131;
t351 = Ifges(4,5) * qJD(3) + qJD(1) * t217 + t183 * t55;
t250 = qJD(1) * qJD(2);
t163 = qJDD(1) * qJ(2) + t250;
t159 = qJDD(1) * t187 + qJDD(2);
t258 = qJD(3) * t180;
t84 = t159 * t184 - t162 * t258;
t256 = qJD(3) * t184;
t85 = t159 * t180 + t162 * t256;
t350 = -t180 * t85 - t184 * t84;
t349 = -m(4) - t354;
t222 = mrSges(4,1) * t184 - mrSges(4,2) * t180;
t286 = Ifges(4,4) * t184;
t348 = qJ(2) * t222 + (-Ifges(4,1) * t180 - t286) * t359;
t347 = mrSges(3,2) - mrSges(2,1) - mrSges(4,3);
t221 = mrSges(4,1) * t180 + mrSges(4,2) * t184;
t223 = pkin(3) * t180 - pkin(7) * t184;
t346 = -m(5) * t223 - t180 * t361 + t184 * t360 + mrSges(2,2) - mrSges(3,3) - t221;
t225 = t138 * t182 - t139 * t178;
t17 = qJD(5) * t225 + t178 * t64 + t182 * t63;
t323 = t17 / 0.2e1;
t71 = t138 * t178 + t139 * t182;
t18 = -qJD(5) * t71 - t178 * t63 + t182 * t64;
t322 = t18 / 0.2e1;
t315 = t63 / 0.2e1;
t314 = t64 / 0.2e1;
t150 = -qJDD(1) * t180 - t184 * t249;
t135 = qJDD(4) - t150;
t129 = qJDD(5) + t135;
t308 = t129 / 0.2e1;
t307 = t135 / 0.2e1;
t157 = t186 * t179;
t158 = t186 * t183;
t82 = t157 * t182 + t158 * t178;
t344 = qJD(5) * t82 + t178 * t357 + t182 * t358;
t83 = t157 * t178 - t158 * t182;
t343 = -qJD(5) * t83 - t178 * t358 + t182 * t357;
t152 = qJ(2) + t223;
t118 = t152 * qJD(1);
t151 = t180 * t162;
t125 = qJD(3) * pkin(7) + t151;
t58 = t118 * t179 + t125 * t183;
t342 = t179 * t58;
t341 = -mrSges(5,1) - t324;
t339 = pkin(4) * t356 - t151;
t107 = t205 * t184;
t143 = t178 * t183 + t179 * t182;
t117 = t143 * qJD(1);
t74 = t330 * t143;
t338 = -qJD(3) * t107 - t180 * t74 - t117;
t105 = t143 * t184;
t337 = qJD(1) * t205 - qJD(3) * t105 + t180 * t189;
t336 = t205 * t180;
t265 = t180 * t187;
t160 = t183 * t265;
t90 = t152 * t179 + t160;
t177 = qJ(4) + qJ(5);
t171 = sin(t177);
t172 = cos(t177);
t220 = -mrSges(5,1) * t183 + mrSges(5,2) * t179;
t335 = m(5) * pkin(3) + mrSges(6,1) * t172 - mrSges(6,2) * t171 - t220 + t361;
t334 = -m(5) * pkin(7) - t360;
t44 = mrSges(5,1) * t135 - mrSges(5,3) * t63;
t45 = -mrSges(5,2) * t135 + mrSges(5,3) * t64;
t333 = -t179 * t44 + t183 * t45;
t253 = qJD(4) * t183;
t66 = -pkin(3) * t150 - pkin(7) * t149 + t163;
t81 = qJDD(3) * pkin(7) + t85;
t21 = t118 * t253 - t125 * t254 + t179 * t66 + t183 * t81;
t22 = -qJD(4) * t58 - t179 * t81 + t183 * t66;
t332 = -t179 * t22 + t183 * t21;
t299 = pkin(4) * t179;
t331 = -t299 + t187;
t49 = pkin(8) * t138 + t58;
t277 = t182 * t49;
t57 = t118 * t183 - t125 * t179;
t48 = -pkin(8) * t139 + t57;
t41 = pkin(4) * t164 + t48;
t14 = t178 * t41 + t277;
t271 = t162 * t184;
t126 = -qJD(3) * pkin(3) - t271;
t79 = -pkin(4) * t138 + t126;
t329 = -mrSges(6,1) * t79 + mrSges(6,3) * t14;
t281 = t178 * t49;
t13 = t182 * t41 - t281;
t328 = mrSges(6,2) * t79 - mrSges(6,3) * t13;
t327 = qJD(1) ^ 2;
t326 = Ifges(6,4) * t323 + Ifges(6,2) * t322 + Ifges(6,6) * t308;
t325 = Ifges(6,1) * t323 + Ifges(6,4) * t322 + Ifges(6,5) * t308;
t321 = Ifges(5,1) * t315 + Ifges(5,4) * t314 + Ifges(5,5) * t307;
t161 = qJD(5) + t164;
t301 = Ifges(6,4) * t71;
t28 = Ifges(6,2) * t225 + Ifges(6,6) * t161 + t301;
t320 = -t28 / 0.2e1;
t319 = t28 / 0.2e1;
t65 = Ifges(6,4) * t225;
t29 = Ifges(6,1) * t71 + Ifges(6,5) * t161 + t65;
t318 = -t29 / 0.2e1;
t317 = t29 / 0.2e1;
t313 = -t225 / 0.2e1;
t312 = t225 / 0.2e1;
t311 = -t71 / 0.2e1;
t310 = t71 / 0.2e1;
t309 = -m(3) - m(4);
t305 = t139 / 0.2e1;
t304 = -t161 / 0.2e1;
t303 = t161 / 0.2e1;
t300 = pkin(4) * t139;
t296 = g(3) * t184;
t268 = t180 * t181;
t100 = t171 * t185 + t172 * t268;
t99 = -t171 * t268 + t172 * t185;
t291 = mrSges(6,1) * t99 - mrSges(6,2) * t100;
t266 = t180 * t185;
t101 = t171 * t266 + t172 * t181;
t102 = -t171 * t181 + t172 * t266;
t290 = mrSges(6,1) * t101 + mrSges(6,2) * t102;
t289 = mrSges(5,3) * t138;
t288 = mrSges(5,3) * t139;
t285 = Ifges(5,4) * t139;
t284 = Ifges(5,4) * t179;
t283 = Ifges(5,4) * t183;
t269 = t179 * t185;
t264 = t181 * t183;
t262 = t183 * t185;
t260 = pkin(1) * t185 + qJ(2) * t181;
t255 = qJD(3) * t187;
t252 = qJD(4) * t184;
t251 = qJDD(1) * mrSges(3,2);
t246 = Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t129;
t245 = Ifges(5,5) * t63 + Ifges(5,6) * t64 + Ifges(5,3) * t135;
t241 = t179 * t265;
t237 = t179 * t258;
t236 = t180 * t255;
t235 = t184 * t255;
t228 = -t179 * t187 + pkin(4);
t227 = -t249 / 0.2e1;
t226 = (t163 + t250) * qJ(2);
t219 = mrSges(5,1) * t179 + mrSges(5,2) * t183;
t218 = -mrSges(6,1) * t171 - mrSges(6,2) * t172;
t216 = Ifges(5,1) * t183 - t284;
t215 = Ifges(5,1) * t179 + t283;
t214 = -t180 * Ifges(4,2) + t286;
t213 = -Ifges(5,2) * t179 + t283;
t212 = Ifges(5,2) * t183 + t284;
t211 = -Ifges(4,5) * t180 - Ifges(4,6) * t184;
t210 = Ifges(5,5) * t183 - Ifges(5,6) * t179;
t209 = Ifges(5,5) * t179 + Ifges(5,6) * t183;
t134 = t183 * t152;
t67 = -pkin(8) * t263 + t180 * t228 + t134;
t78 = -pkin(8) * t270 + t90;
t32 = -t178 * t78 + t182 * t67;
t33 = t178 * t67 + t182 * t78;
t208 = t183 * t57 + t342;
t87 = -mrSges(5,2) * t164 + t289;
t88 = mrSges(5,1) * t164 - t288;
t207 = -t179 * t87 - t183 * t88;
t10 = pkin(8) * t64 + t21;
t9 = pkin(4) * t135 - pkin(8) * t63 + t22;
t2 = qJD(5) * t13 + t10 * t182 + t178 * t9;
t3 = -qJD(5) * t14 - t10 * t178 + t182 * t9;
t204 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t246;
t121 = t179 * t266 + t264;
t119 = -t179 * t268 + t262;
t80 = -qJDD(3) * pkin(3) - t84;
t202 = t180 * (-Ifges(4,2) * t184 - t287);
t197 = t179 * t252 + t180 * t257;
t196 = -t183 * t252 + t237;
t136 = qJD(3) * t224 + qJD(2);
t46 = -qJD(4) * t241 + t136 * t179 + t152 * t253 + t183 * t235;
t193 = Ifges(5,5) * t184 - t180 * t216;
t192 = Ifges(5,6) * t184 - t180 * t213;
t191 = Ifges(5,3) * t184 - t180 * t210;
t190 = -qJD(4) * t208 + t332;
t167 = -pkin(1) * qJDD(1) + qJDD(2);
t154 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t170;
t144 = t221 * qJD(1);
t137 = t331 * t184;
t130 = t219 * t184;
t122 = -t179 * t181 + t180 * t262;
t120 = t180 * t264 + t269;
t114 = Ifges(4,6) * qJD(3) + qJD(1) * t214;
t111 = t183 * t136;
t109 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t150;
t104 = t143 * t180;
t96 = qJD(1) * t336;
t95 = t180 * t117;
t89 = t134 - t241;
t86 = -pkin(4) * t196 + t236;
t54 = Ifges(5,2) * t138 + Ifges(5,6) * t164 + t285;
t53 = Ifges(5,5) * t139 + Ifges(5,6) * t138 + Ifges(5,3) * t164;
t51 = mrSges(6,1) * t161 - mrSges(6,3) * t71;
t50 = -mrSges(6,2) * t161 + mrSges(6,3) * t225;
t47 = -qJD(4) * t90 - t179 * t235 + t111;
t40 = t143 * t258 + t184 * t189;
t38 = qJD(3) * t336 - t184 * t74;
t36 = -pkin(4) * t64 + t80;
t35 = pkin(8) * t196 + t46;
t34 = -mrSges(6,1) * t225 + mrSges(6,2) * t71;
t31 = t111 + (-t160 + (pkin(8) * t184 - t152) * t179) * qJD(4) + (t184 * t228 + t247) * qJD(3);
t27 = Ifges(6,5) * t71 + Ifges(6,6) * t225 + Ifges(6,3) * t161;
t23 = Ifges(5,4) * t63 + Ifges(5,2) * t64 + Ifges(5,6) * t135;
t20 = t182 * t48 - t281;
t19 = -t178 * t48 - t277;
t12 = -mrSges(6,2) * t129 + mrSges(6,3) * t18;
t11 = mrSges(6,1) * t129 - mrSges(6,3) * t17;
t8 = -qJD(5) * t33 - t178 * t35 + t182 * t31;
t7 = qJD(5) * t32 + t178 * t31 + t182 * t35;
t6 = -mrSges(6,1) * t18 + mrSges(6,2) * t17;
t1 = [(Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (t221 + 0.2e1 * mrSges(3,3)) * t163 + (Ifges(4,1) * t149 + Ifges(4,4) * t150) * t359 + t154 * t235 + (-Ifges(6,4) * t107 - Ifges(6,2) * t105 + Ifges(6,6) * t180) * t322 + (-Ifges(6,1) * t107 - Ifges(6,4) * t105 + Ifges(6,5) * t180) * t323 + (-Ifges(6,5) * t107 - Ifges(6,6) * t105 + Ifges(6,3) * t180) * t308 + t36 * (mrSges(6,1) * t105 - mrSges(6,2) * t107) + t3 * (mrSges(6,1) * t180 + mrSges(6,3) * t107) - t351 * t258 / 0.2e1 + t348 * t249 + t350 * mrSges(4,3) - (t179 * t55 + t183 * t54) * t252 / 0.2e1 + (t53 + t27) * t256 / 0.2e1 + (t246 + t245) * t180 / 0.2e1 + m(6) * (t13 * t8 - t137 * t36 + t14 * t7 + t2 * t33 + t3 * t32 + t79 * t86) + t202 * t227 + t32 * t11 + t33 * t12 + t126 * (-mrSges(5,1) * t196 - mrSges(5,2) * t197) + (-t122 * mrSges(5,1) - t102 * mrSges(6,1) + t121 * mrSges(5,2) + t101 * mrSges(6,2) + (m(3) * pkin(1) - t331 * m(6) + (-m(4) - m(5)) * t187 - t347) * t181 + ((-m(3) + t349) * qJ(2) + t346) * t185) * g(1) + t353 * t184 * t187 + m(4) * (-t187 * t350 + t226) - t180 * (Ifges(4,4) * t149 + Ifges(4,2) * t150) / 0.2e1 - t107 * t325 - t105 * t326 + (Ifges(6,4) * t38 + Ifges(6,2) * t40 + Ifges(6,6) * t256) * t312 + (Ifges(5,6) * t180 + t184 * t213) * t314 + (Ifges(5,5) * t180 + t184 * t216) * t315 + t38 * t317 + t40 * t319 + t263 * t321 + (Ifges(6,5) * t38 + Ifges(6,6) * t40 + Ifges(6,3) * t256) * t303 + (qJD(3) * t193 - t215 * t252) * t305 + (Ifges(5,3) * t180 + t184 * t210) * t307 + (Ifges(6,1) * t38 + Ifges(6,4) * t40 + Ifges(6,5) * t256) * t310 + t2 * (-mrSges(6,2) * t180 - mrSges(6,3) * t105) + qJD(3) ^ 2 * t211 / 0.2e1 + t150 * t214 / 0.2e1 + t149 * t217 / 0.2e1 + (-t269 * t324 - m(3) * t260 - t120 * mrSges(5,1) - t100 * mrSges(6,1) - t119 * mrSges(5,2) - t99 * mrSges(6,2) + t349 * (pkin(6) * t185 + t260) + t347 * t185 + t346 * t181) * g(2) + t109 * t265 - t23 * t270 / 0.2e1 + t21 * (-mrSges(5,2) * t180 - mrSges(5,3) * t270) + t22 * (mrSges(5,1) * t180 - mrSges(5,3) * t263) + t352 * t236 + t7 * t50 + t8 * t51 + t79 * (-mrSges(6,1) * t40 + mrSges(6,2) * t38) + t86 * t34 + t46 * t87 + t47 * t88 + t89 * t44 + t90 * t45 + t80 * t130 - t137 * t6 + qJD(2) * t144 + qJ(2) * (-mrSges(4,1) * t150 + mrSges(4,2) * t149) + t167 * mrSges(3,2) + qJDD(3) * (Ifges(4,5) * t184 - Ifges(4,6) * t180) + m(3) * (-pkin(1) * t167 + t226) + t54 * t237 / 0.2e1 - pkin(1) * t251 + t138 * (qJD(3) * t192 - t212 * t252) / 0.2e1 + t164 * (qJD(3) * t191 - t209 * t252) / 0.2e1 - t114 * t256 / 0.2e1 + t14 * (-mrSges(6,2) * t256 + mrSges(6,3) * t40) + t58 * (-mrSges(5,2) * t256 + mrSges(5,3) * t196) + t13 * (mrSges(6,1) * t256 - mrSges(6,3) * t38) + t57 * (mrSges(5,1) * t256 + mrSges(5,3) * t197) + m(5) * (t21 * t90 + t22 * t89 + t46 * t58 + t47 * t57 + (t126 * t258 - t184 * t80) * t187); t251 - t104 * t11 - t336 * t12 + t337 * t51 + t338 * t50 + (qJ(2) * t309 - mrSges(3,3)) * t327 + (-t144 + t207) * qJD(1) + (-t6 + (-t179 * t88 + t183 * t87 + t154) * qJD(3) + t353) * t184 + (t109 + t207 * qJD(4) + (t34 + t352) * qJD(3) + t333) * t180 + m(3) * t167 - m(4) * t350 - t355 * (-t309 + t354) + (-t104 * t3 + t13 * t337 + t14 * t338 - t184 * t36 - t2 * t336 + t258 * t79) * m(6) + ((-t80 + (-t179 * t57 + t183 * t58) * qJD(3)) * t184 + (qJD(3) * t126 + t190) * t180 - t208 * qJD(1)) * m(5); t351 * t170 / 0.2e1 - (Ifges(6,1) * t310 + Ifges(6,4) * t312 + Ifges(6,5) * t303 + t317 + t328) * t189 + (t202 / 0.2e1 - t348) * t327 + (t13 * t96 - t14 * t95 - t143 * t3 - t2 * t205) * mrSges(6,3) + (Ifges(6,4) * t143 - Ifges(6,2) * t205) * t322 + (Ifges(6,1) * t143 - Ifges(6,4) * t205) * t323 + (Ifges(6,5) * t143 - Ifges(6,6) * t205) * t308 + t36 * (mrSges(6,1) * t205 + mrSges(6,2) * t143) - t205 * t326 + (t180 * t335 + t184 * t334 + t221) * g(3) + t211 * t227 - pkin(3) * t30 + (Ifges(6,4) * t96 + Ifges(6,2) * t95) * t313 - (Ifges(6,4) * t310 + Ifges(6,2) * t312 + Ifges(6,6) * t303 + t319 + t329) * t74 + (-mrSges(5,1) * t57 + mrSges(5,2) * t58 - t53 / 0.2e1 - t27 / 0.2e1 + Ifges(6,5) * t311 + Ifges(6,6) * t313 + Ifges(6,3) * t304 + t14 * mrSges(6,2) - t13 * mrSges(6,1) + t114 / 0.2e1) * t259 + (Ifges(6,1) * t96 + Ifges(6,4) * t95) * t311 + t164 * t219 * t126 + t343 * t51 + (t13 * t343 + t14 * t344 - t166 * t36 + t2 * t83 + t3 * t82 + t339 * t79) * m(6) + t344 * t50 - t352 * t151 + (-pkin(3) * t80 - t126 * t151 - t57 * t75 - t58 * t76) * m(5) + t143 * t325 + t212 * t314 + t215 * t315 + t96 * t318 + t95 * t320 + t179 * t321 + t209 * t307 + (-t253 * t57 - t254 * t58 + (-t180 * t342 - t267 * t57) * qJD(1) + t332) * mrSges(5,3) + t355 * (t180 * t334 - t184 * t335 - t222) - t356 * t54 / 0.2e1 + (t138 * t213 + t139 * t216 + t164 * t210) * qJD(4) / 0.2e1 - (t138 * t192 + t139 * t193 + t164 * t191) * qJD(1) / 0.2e1 + t80 * t220 + (m(5) * t190 - t253 * t88 - t254 * t87 + t333) * pkin(7) + (Ifges(6,5) * t96 + Ifges(6,6) * t95) * t304 - t154 * t271 + t82 * t11 + t83 * t12 + t84 * mrSges(4,1) - t85 * mrSges(4,2) + t339 * t34 - t76 * t87 - t75 * t88 - t79 * (-mrSges(6,1) * t95 + mrSges(6,2) * t96) + Ifges(4,5) * t149 + Ifges(4,6) * t150 - t166 * t6 + t183 * t23 / 0.2e1 + Ifges(4,3) * qJDD(3) + t55 * t253 / 0.2e1; t22 * mrSges(5,1) - t21 * mrSges(5,2) - t139 * (Ifges(5,1) * t138 - t285) / 0.2e1 + t245 + (t178 * t2 + t182 * t3 + (-t13 * t178 + t14 * t182) * qJD(5)) * t324 + t54 * t305 - (-m(6) * t299 + t218) * t296 - t34 * t300 - m(6) * (t13 * t19 + t14 * t20 + t300 * t79) + t204 - t20 * t50 - t19 * t51 + g(3) * t130 - t126 * (mrSges(5,1) * t139 + mrSges(5,2) * t138) - t164 * (Ifges(5,5) * t138 - Ifges(5,6) * t139) / 0.2e1 + (Ifges(6,1) * t311 + Ifges(6,4) * t313 + Ifges(6,5) * t304 + t318 - t328) * t225 - (Ifges(6,4) * t311 + Ifges(6,2) * t313 + Ifges(6,6) * t304 + t320 - t329) * t71 + (t288 + t88) * t58 + (t289 - t87) * t57 - (-Ifges(5,2) * t139 + t131 + t55) * t138 / 0.2e1 + (-mrSges(5,2) * t122 + t121 * t341 - t290) * g(2) + (mrSges(5,2) * t120 + t119 * t341 - t291) * g(1) + ((-t178 * t51 + t182 * t50) * qJD(5) + t182 * t11 + t178 * t12) * pkin(4); -t79 * (mrSges(6,1) * t71 + mrSges(6,2) * t225) + (Ifges(6,1) * t225 - t301) * t311 + t28 * t310 + (Ifges(6,5) * t225 - Ifges(6,6) * t71) * t304 - t13 * t50 + t14 * t51 - g(1) * t291 - g(2) * t290 - t218 * t296 + (t13 * t225 + t14 * t71) * mrSges(6,3) + t204 + (-Ifges(6,2) * t71 + t29 + t65) * t313;];
tau = t1;
