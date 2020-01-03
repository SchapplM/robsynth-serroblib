% Calculate vector of inverse dynamics joint torques for
% S5RRRPR3
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:51
% EndTime: 2020-01-03 12:09:02
% DurationCPUTime: 4.94s
% Computational Cost: add. (6393->452), mult. (9682->602), div. (0->0), fcn. (6446->16), ass. (0->234)
t267 = sin(qJ(3));
t351 = t267 / 0.2e1;
t271 = cos(qJ(3));
t248 = t271 * qJD(4);
t265 = -qJ(4) - pkin(7);
t294 = qJD(3) * t265;
t175 = t267 * t294 + t248;
t176 = -qJD(4) * t267 + t271 * t294;
t263 = sin(pkin(9));
t264 = cos(pkin(9));
t190 = t263 * t271 + t264 * t267;
t272 = cos(qJ(2));
t330 = qJD(1) * pkin(1);
t301 = t272 * t330;
t364 = -t175 * t263 + t264 * t176 + t190 * t301;
t189 = -t263 * t267 + t264 * t271;
t363 = t264 * t175 + t263 * t176 - t189 * t301;
t214 = -mrSges(4,1) * t271 + t267 * mrSges(4,2);
t385 = -mrSges(3,1) + t214;
t179 = t189 * qJD(3);
t342 = pkin(8) * t179;
t384 = -t342 + t364;
t178 = t190 * qJD(3);
t170 = t178 * pkin(8);
t383 = t170 - t363;
t259 = qJ(3) + pkin(9);
t245 = sin(t259);
t246 = cos(t259);
t247 = qJ(5) + t259;
t231 = sin(t247);
t328 = t231 * mrSges(6,2);
t382 = -t246 * mrSges(5,1) + t245 * mrSges(5,2) + t328;
t266 = sin(qJ(5));
t270 = cos(qJ(5));
t258 = qJD(1) + qJD(2);
t160 = t190 * t258;
t343 = pkin(8) * t160;
t268 = sin(qJ(2));
t302 = t268 * t330;
t203 = pkin(7) * t258 + t302;
t293 = qJ(4) * t258 + t203;
t145 = t293 * t271;
t129 = t263 * t145;
t144 = t293 * t267;
t329 = qJD(3) * pkin(3);
t135 = -t144 + t329;
t75 = t264 * t135 - t129;
t50 = qJD(3) * pkin(4) - t343 + t75;
t159 = t189 * t258;
t344 = pkin(8) * t159;
t312 = t264 * t145;
t76 = t263 * t135 + t312;
t52 = t76 + t344;
t22 = -t266 * t52 + t270 * t50;
t23 = t266 * t50 + t270 * t52;
t262 = qJ(1) + qJ(2);
t250 = cos(t262);
t254 = qJDD(3) + qJDD(5);
t257 = qJD(3) + qJD(5);
t292 = t270 * t159 - t160 * t266;
t255 = qJDD(1) + qJDD(2);
t304 = qJD(3) * t267;
t180 = t255 * t271 - t258 * t304;
t303 = qJD(3) * t271;
t181 = t255 * t267 + t258 * t303;
t106 = t180 * t263 + t181 * t264;
t298 = qJD(2) * t330;
t323 = qJDD(1) * pkin(1);
t188 = t268 * t323 + t272 * t298;
t166 = pkin(7) * t255 + t188;
t297 = t203 * t303;
t65 = -t297 + qJDD(3) * pkin(3) - qJ(4) * t181 + (-qJD(4) * t258 - t166) * t267;
t108 = t271 * t166 - t203 * t304;
t70 = qJ(4) * t180 + t248 * t258 + t108;
t31 = -t263 * t70 + t264 * t65;
t18 = qJDD(3) * pkin(4) - pkin(8) * t106 + t31;
t105 = t180 * t264 - t181 * t263;
t32 = t263 * t65 + t264 * t70;
t19 = pkin(8) * t105 + t32;
t3 = qJD(5) * t22 + t18 * t266 + t19 * t270;
t232 = cos(t247);
t317 = t232 * t250;
t336 = mrSges(6,1) * t231;
t92 = t159 * t266 + t160 * t270;
t349 = Ifges(6,4) * t92;
t36 = qJD(5) * t292 + t105 * t266 + t106 * t270;
t37 = -qJD(5) * t92 + t105 * t270 - t106 * t266;
t4 = -qJD(5) * t23 + t18 * t270 - t19 * t266;
t86 = Ifges(6,4) * t292;
t41 = Ifges(6,1) * t92 + Ifges(6,5) * t257 + t86;
t345 = pkin(3) * t271;
t238 = pkin(2) + t345;
t156 = -t238 * t258 + qJD(4) - t301;
t99 = -pkin(4) * t159 + t156;
t381 = t4 * mrSges(6,1) + Ifges(6,3) * t254 + Ifges(6,6) * t37 + Ifges(6,5) * t36 - g(3) * (mrSges(6,2) * t317 + t250 * t336) - t3 * mrSges(6,2) - (Ifges(6,5) * t292 - Ifges(6,6) * t92) * t257 / 0.2e1 + (t22 * t292 + t23 * t92) * mrSges(6,3) - (-Ifges(6,2) * t92 + t41 + t86) * t292 / 0.2e1 - t99 * (mrSges(6,1) * t92 + mrSges(6,2) * t292) - (Ifges(6,1) * t292 - t349) * t92 / 0.2e1;
t338 = qJD(3) / 0.2e1;
t233 = pkin(3) * t264 + pkin(4);
t347 = pkin(3) * t263;
t174 = t233 * t266 + t270 * t347;
t80 = t144 * t263 - t312;
t54 = t80 - t344;
t81 = -t264 * t144 - t129;
t55 = t81 - t343;
t377 = -t174 * qJD(5) + t266 * t55 - t270 * t54;
t173 = t233 * t270 - t266 * t347;
t376 = t173 * qJD(5) - t266 * t54 - t270 * t55;
t315 = t258 * t267;
t201 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t315;
t314 = t258 * t271;
t202 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t314;
t375 = (t258 * mrSges(3,2) + t267 * t201 - t271 * t202) * t272;
t374 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t373 = t382 + t385;
t40 = Ifges(6,2) * t292 + Ifges(6,6) * t257 + t349;
t371 = t40 / 0.2e1;
t213 = t265 * t267;
t251 = t271 * qJ(4);
t215 = pkin(7) * t271 + t251;
t133 = t264 * t213 - t215 * t263;
t341 = pkin(8) * t190;
t103 = t133 - t341;
t134 = t263 * t213 + t264 * t215;
t184 = t189 * pkin(8);
t104 = t184 + t134;
t45 = t103 * t270 - t104 * t266;
t370 = qJD(5) * t45 + t266 * t384 - t270 * t383;
t46 = t103 * t266 + t104 * t270;
t369 = -qJD(5) * t46 + t266 * t383 + t270 * t384;
t362 = t385 * t258;
t335 = mrSges(4,2) * t271;
t346 = pkin(3) * t267;
t358 = m(5) * pkin(3);
t361 = -t267 * (mrSges(4,1) + t358) + m(6) * (-pkin(4) * t245 - t346) - mrSges(5,1) * t245 - mrSges(5,2) * t246 - t335;
t204 = -pkin(2) * t258 - t301;
t359 = t204 * t268 + (t267 ^ 2 + t271 ^ 2) * t203 * t272;
t355 = t92 / 0.2e1;
t353 = t160 / 0.2e1;
t348 = pkin(1) * t272;
t249 = sin(t262);
t340 = g(2) * t249;
t339 = -qJD(3) / 0.2e1;
t334 = Ifges(4,4) * t267;
t333 = Ifges(4,4) * t271;
t332 = Ifges(5,4) * t160;
t331 = pkin(1) * qJD(2);
t221 = t232 * mrSges(6,1);
t325 = t271 * Ifges(4,2);
t322 = t108 * t271;
t109 = -t166 * t267 - t297;
t321 = t109 * t267;
t237 = pkin(1) * t268 + pkin(7);
t316 = t237 * t271;
t308 = -qJ(4) - t237;
t289 = qJD(3) * t308;
t300 = t272 * t331;
t121 = t267 * t289 + t271 * t300 + t248;
t122 = (-qJD(4) - t300) * t267 + t271 * t289;
t67 = t264 * t121 + t263 * t122;
t185 = t308 * t267;
t186 = t251 + t316;
t112 = t263 * t185 + t264 * t186;
t295 = pkin(4) * t246 + t345;
t199 = pkin(2) + t295;
t256 = -pkin(8) + t265;
t307 = t249 * t199 + t250 * t256;
t306 = t249 * t238 + t250 * t265;
t305 = t250 * pkin(2) + t249 * pkin(7);
t241 = pkin(3) * t304;
t9 = -t37 * mrSges(6,1) + t36 * mrSges(6,2);
t296 = t303 / 0.2e1;
t142 = pkin(4) * t178 + t241;
t49 = -t105 * mrSges(5,1) + t106 * mrSges(5,2);
t66 = -t121 * t263 + t264 * t122;
t111 = t264 * t185 - t186 * t263;
t291 = t250 * t199 - t249 * t256;
t290 = t250 * t238 - t249 * t265;
t187 = -t268 * t298 + t272 * t323;
t288 = g(2) * t250 + g(3) * t249;
t287 = -mrSges(6,2) * t232 - t336;
t285 = t325 + t334;
t284 = Ifges(4,5) * t271 - Ifges(4,6) * t267;
t84 = t111 - t341;
t85 = t184 + t112;
t38 = -t266 * t85 + t270 * t84;
t39 = t266 * t84 + t270 * t85;
t115 = t189 * t270 - t190 * t266;
t116 = t189 * t266 + t190 * t270;
t283 = t201 * t271 + t202 * t267;
t152 = -pkin(4) * t189 - t238;
t282 = t204 * (mrSges(4,1) * t267 + t335);
t281 = t267 * (Ifges(4,1) * t271 - t334);
t165 = -pkin(2) * t255 - t187;
t110 = -pkin(3) * t180 + qJDD(4) + t165;
t277 = -qJD(3) * t283 + t271 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t180) - t267 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t181);
t276 = -mrSges(6,1) * t317 + t374 * t249 + t373 * t250;
t275 = (m(4) * pkin(7) - t374) * t250 + (-t221 + t373) * t249;
t163 = Ifges(4,6) * qJD(3) + t258 * t285;
t220 = Ifges(4,4) * t314;
t164 = Ifges(4,1) * t315 + Ifges(4,5) * qJD(3) + t220;
t51 = -pkin(4) * t105 + t110;
t57 = qJD(5) * t115 - t178 * t266 + t179 * t270;
t58 = -qJD(5) * t116 - t178 * t270 - t179 * t266;
t87 = Ifges(5,2) * t159 + Ifges(5,6) * qJD(3) + t332;
t149 = Ifges(5,4) * t159;
t88 = Ifges(5,1) * t160 + Ifges(5,5) * qJD(3) + t149;
t274 = (t284 * t338 + t282) * qJD(3) + t180 * t285 / 0.2e1 + (Ifges(4,1) * t181 + Ifges(4,4) * t180) * t351 + t181 * (t267 * Ifges(4,1) + t333) / 0.2e1 + (-mrSges(6,1) * t51 + mrSges(6,3) * t3 + Ifges(6,4) * t36 + Ifges(6,2) * t37 + Ifges(6,6) * t254) * t115 + t58 * t371 + (Ifges(5,1) * t179 - Ifges(5,4) * t178) * t353 + (Ifges(5,5) * t179 - Ifges(5,6) * t178) * t338 + t156 * (mrSges(5,1) * t178 + mrSges(5,2) * t179) + t159 * (Ifges(5,4) * t179 - Ifges(5,2) * t178) / 0.2e1 + t292 * (Ifges(6,4) * t57 + Ifges(6,2) * t58) / 0.2e1 + t164 * t296 + (Ifges(6,1) * t57 + Ifges(6,4) * t58) * t355 + mrSges(4,3) * t322 + ((-Ifges(4,2) * t267 + t333) * t296 + t281 * t338) * t258 - t163 * t304 / 0.2e1 + t271 * (Ifges(4,4) * t181 + Ifges(4,2) * t180) / 0.2e1 + (0.2e1 * Ifges(4,5) * t351 + Ifges(5,5) * t190 + Ifges(4,6) * t271 + Ifges(5,6) * t189) * qJDD(3) + (mrSges(6,2) * t51 - mrSges(6,3) * t4 + Ifges(6,1) * t36 + Ifges(6,4) * t37 + Ifges(6,5) * t254) * t116 + t57 * t41 / 0.2e1 + t99 * (-mrSges(6,1) * t58 + mrSges(6,2) * t57) - t178 * t87 / 0.2e1 + t179 * t88 / 0.2e1 + t187 * mrSges(3,1) - t188 * mrSges(3,2) + t110 * (-mrSges(5,1) * t189 + mrSges(5,2) * t190) + t106 * (Ifges(5,1) * t190 + Ifges(5,4) * t189) + t105 * (Ifges(5,4) * t190 + Ifges(5,2) * t189) + t165 * t214 + (-t178 * t76 - t179 * t75 + t189 * t32 - t190 * t31) * mrSges(5,3) + Ifges(3,3) * t255 + t257 * (Ifges(6,5) * t57 + Ifges(6,6) * t58) / 0.2e1 + (-t22 * t57 + t23 * t58) * mrSges(6,3);
t273 = cos(qJ(1));
t269 = sin(qJ(1));
t253 = t273 * pkin(1);
t252 = t269 * pkin(1);
t242 = t268 * t331;
t239 = -pkin(2) - t348;
t235 = t249 * pkin(2);
t210 = -t238 - t348;
t200 = t242 + t241;
t140 = t152 - t348;
t137 = qJD(3) * mrSges(5,1) - t160 * mrSges(5,3);
t136 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t159;
t125 = t142 + t242;
t120 = pkin(3) * t315 + pkin(4) * t160;
t113 = -mrSges(4,1) * t180 + mrSges(4,2) * t181;
t95 = -mrSges(5,1) * t159 + mrSges(5,2) * t160;
t94 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t106;
t93 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t105;
t79 = mrSges(6,1) * t257 - mrSges(6,3) * t92;
t78 = -mrSges(6,2) * t257 + mrSges(6,3) * t292;
t48 = -t170 + t67;
t47 = t66 - t342;
t44 = -mrSges(6,1) * t292 + mrSges(6,2) * t92;
t30 = -mrSges(6,2) * t254 + mrSges(6,3) * t37;
t29 = mrSges(6,1) * t254 - mrSges(6,3) * t36;
t8 = -qJD(5) * t39 - t266 * t48 + t270 * t47;
t7 = qJD(5) * t38 + t266 * t47 + t270 * t48;
t1 = [m(4) * (t108 * t316 + t165 * t239 - t237 * t321) + t277 * t237 + ((mrSges(3,1) * t272 - mrSges(3,2) * t268) * t255 + (-g(2) * t273 - g(3) * t269 + t187 * t272 + t188 * t268) * m(3) + (m(4) * t359 + t268 * t362 - t375) * qJD(2)) * pkin(1) + m(6) * (t125 * t99 + t140 * t51 + t22 * t8 + t23 * t7 + t3 * t39 + t38 * t4) + m(5) * (t110 * t210 + t111 * t31 + t112 * t32 + t156 * t200 + t66 * t75 + t67 * t76) + t274 + (-m(6) * (t252 + t307) - m(5) * (t252 + t306) - m(4) * (t235 + t252) - t269 * mrSges(2,1) - mrSges(2,2) * t273 + t275) * g(3) + t38 * t29 + t39 * t30 - mrSges(4,3) * t321 + t7 * t78 + t8 * t79 + (-m(6) * (t253 + t291) - m(5) * (t253 + t290) - m(4) * (t253 + t305) - mrSges(2,1) * t273 + t269 * mrSges(2,2) + t276) * g(2) + t111 * t94 + t112 * t93 + t125 * t44 + t67 * t136 + t66 * t137 + t140 * t9 + t200 * t95 + t210 * t49 + t239 * t113 + Ifges(2,3) * qJDD(1); t276 * g(2) + t277 * pkin(7) + t275 * g(3) + t274 + t369 * t79 + t370 * t78 + (t375 + (-t362 - t44 - t95) * t268) * t330 + t364 * t137 + t363 * t136 + (-t109 * mrSges(4,3) + t329 * t95) * t267 + t45 * t29 + t46 * t30 - pkin(2) * t113 + t133 * t94 + t134 * t93 + t142 * t44 + t152 * t9 - t238 * t49 + (-g(2) * t291 - g(3) * t307 + t152 * t51 + t3 * t46 + t4 * t45 + (t142 - t302) * t99 + t370 * t23 + t369 * t22) * m(6) + (-g(2) * t290 - g(3) * t306 - t110 * t238 + t133 * t31 + t134 * t32 + t363 * t76 + t364 * t75 + (t241 - t302) * t156) * m(5) + (-t305 * g(2) - t235 * g(3) - pkin(2) * t165 + (-t321 + t322) * pkin(7) - t359 * t330) * m(4); t92 * t371 - (-Ifges(5,2) * t160 + t149 + t88) * t159 / 0.2e1 + (t159 * t75 + t160 * t76) * mrSges(5,3) + (Ifges(5,5) * t159 - Ifges(5,6) * t160) * t339 - t156 * (mrSges(5,1) * t160 + mrSges(5,2) * t159) - t160 * (Ifges(5,1) * t159 - t332) / 0.2e1 + t283 * t203 + (-m(5) * t345 + t214 - t221 + t382) * g(1) - m(5) * (t75 * t80 + t76 * t81) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t163 * t351 - t282 + t284 * t339 + (t325 * t351 - t281 / 0.2e1) * t258 + (-m(5) * t156 - t95) * t346 - (t164 + t220) * t271 / 0.2e1) * t258 + t87 * t353 + (t263 * t32 + t264 * t31) * t358 + t381 + (-t287 - t361) * t340 + t361 * g(3) * t250 + t31 * mrSges(5,1) - t32 * mrSges(5,2) + Ifges(5,6) * t105 + Ifges(5,5) * t106 - t108 * mrSges(4,2) + t109 * mrSges(4,1) - t120 * t44 - t81 * t136 - t80 * t137 + t173 * t29 + t174 * t30 + Ifges(4,6) * t180 + Ifges(4,5) * t181 + (t263 * t93 + t264 * t94) * pkin(3) + t376 * t78 + t377 * t79 + (-t295 * g(1) - t120 * t99 + t173 * t4 + t174 * t3 + t377 * t22 + t376 * t23) * m(6); -t159 * t136 + t160 * t137 - t292 * t78 + t92 * t79 + t49 + t9 + (t22 * t92 - t23 * t292 + t288 + t51) * m(6) + (-t159 * t76 + t160 * t75 + t110 + t288) * m(5); t40 * t355 - t22 * t78 + t23 * t79 - g(1) * (t221 - t328) - t287 * t340 + t381;];
tau = t1;
