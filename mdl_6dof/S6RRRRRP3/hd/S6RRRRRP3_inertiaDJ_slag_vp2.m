% Calculate time derivative of joint inertia matrix for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:05:54
% EndTime: 2019-03-10 01:06:08
% DurationCPUTime: 6.06s
% Computational Cost: add. (10333->515), mult. (23096->695), div. (0->0), fcn. (21695->8), ass. (0->215)
t253 = sin(qJ(3));
t254 = sin(qJ(2));
t257 = cos(qJ(3));
t258 = cos(qJ(2));
t219 = t253 * t258 + t257 * t254;
t356 = qJD(2) + qJD(3);
t170 = t356 * t219;
t374 = Ifges(6,6) + Ifges(7,6);
t382 = t374 * t170;
t375 = Ifges(6,5) + Ifges(7,5);
t381 = t375 * t170;
t380 = Ifges(6,1) + Ifges(7,1);
t379 = Ifges(6,4) + Ifges(7,4);
t378 = Ifges(6,2) + Ifges(7,2);
t342 = -pkin(8) - pkin(7);
t236 = t342 * t254;
t238 = t342 * t258;
t189 = t236 * t253 - t238 * t257;
t288 = qJD(2) * t342;
t230 = t254 * t288;
t276 = t258 * t288;
t121 = t189 * qJD(3) + t230 * t253 - t257 * t276;
t252 = sin(qJ(4));
t299 = qJD(4) * t252;
t286 = t219 * t299;
t217 = t253 * t254 - t257 * t258;
t169 = t356 * t217;
t256 = cos(qJ(4));
t311 = t169 * t256;
t264 = t286 + t311;
t298 = qJD(4) * t256;
t312 = t169 * t252;
t265 = t219 * t298 - t312;
t66 = t265 * mrSges(5,1) - t264 * mrSges(5,2);
t377 = m(5) * t121 + t66;
t251 = sin(qJ(5));
t255 = cos(qJ(5));
t268 = t251 * t252 - t255 * t256;
t355 = qJD(4) + qJD(5);
t167 = t355 * t268;
t218 = t251 * t256 + t252 * t255;
t168 = t355 * t218;
t277 = -t167 * t375 - t168 * t374;
t376 = Ifges(5,5) * t298 + t277;
t357 = (t252 ^ 2 + t256 ^ 2) * t257;
t373 = Ifges(6,3) + Ifges(7,3);
t48 = -t168 * t219 + t268 * t169;
t143 = t268 * t219;
t49 = t143 * t355 + t218 * t169;
t372 = t378 * t49 + t379 * t48 + t382;
t371 = t379 * t49 + t380 * t48 + t381;
t142 = t218 * t219;
t370 = -t142 * t378 - t143 * t379 + t217 * t374;
t369 = -t379 * t142 - t143 * t380 + t375 * t217;
t368 = -t167 * t379 - t168 * t378;
t367 = -t167 * t380 - t379 * t168;
t366 = t218 * t379 - t268 * t378;
t365 = t218 * t380 - t379 * t268;
t358 = t257 * t236 + t238 * t253;
t120 = qJD(3) * t358 + t257 * t230 + t253 * t276;
t244 = -pkin(2) * t258 - pkin(1);
t153 = t217 * pkin(3) - t219 * pkin(9) + t244;
t99 = pkin(2) * qJD(2) * t254 + pkin(3) * t170 + pkin(9) * t169;
t31 = t256 * t120 + t153 * t298 - t189 * t299 + t252 * t99;
t179 = t256 * t189;
t102 = t252 * t153 + t179;
t280 = -t120 * t252 + t256 * t99;
t32 = -t102 * qJD(4) + t280;
t364 = -t32 * t252 + t256 * t31;
t362 = pkin(4) * t255;
t361 = -mrSges(6,1) - mrSges(7,1);
t360 = pkin(4) * qJD(5);
t359 = -Ifges(5,5) * t311 + Ifges(5,3) * t170;
t240 = pkin(2) * t253 + pkin(9);
t325 = -pkin(10) - t240;
t213 = t325 * t252;
t248 = t256 * pkin(10);
t214 = t240 * t256 + t248;
t147 = t251 * t213 + t255 * t214;
t341 = -pkin(10) - pkin(9);
t235 = t341 * t252;
t237 = pkin(9) * t256 + t248;
t188 = t251 * t235 + t255 * t237;
t302 = -t168 * qJ(6) - qJD(6) * t268;
t281 = qJD(4) * t325;
t321 = pkin(2) * qJD(3);
t291 = t257 * t321;
t190 = t252 * t281 + t256 * t291;
t191 = -t252 * t291 + t256 * t281;
t296 = qJD(5) * t255;
t297 = qJD(5) * t251;
t78 = t255 * t190 + t251 * t191 + t213 * t296 - t214 * t297;
t54 = t78 + t302;
t267 = qJ(6) * t167 - qJD(6) * t218;
t79 = -qJD(5) * t147 - t190 * t251 + t255 * t191;
t55 = t267 + t79;
t354 = t79 * mrSges(6,1) + t55 * mrSges(7,1) - t78 * mrSges(6,2) - t54 * mrSges(7,2);
t287 = qJD(4) * t341;
t227 = t252 * t287;
t228 = t256 * t287;
t118 = t255 * t227 + t251 * t228 + t235 * t296 - t237 * t297;
t119 = -qJD(5) * t188 - t227 * t251 + t255 * t228;
t67 = t118 + t302;
t68 = t119 + t267;
t353 = t119 * mrSges(6,1) + t68 * mrSges(7,1) - t118 * mrSges(6,2) - t67 * mrSges(7,2);
t352 = 0.2e1 * m(5);
t351 = 2 * m(6);
t350 = 2 * m(7);
t103 = t168 * mrSges(7,1) - t167 * mrSges(7,2);
t349 = 0.2e1 * t103;
t104 = mrSges(6,1) * t168 - mrSges(6,2) * t167;
t348 = 0.2e1 * t104;
t347 = 0.2e1 * t121;
t173 = mrSges(7,1) * t268 + mrSges(7,2) * t218;
t346 = 0.2e1 * t173;
t345 = 0.2e1 * t244;
t330 = pkin(2) * t257;
t174 = mrSges(6,1) * t268 + mrSges(6,2) * t218;
t329 = pkin(4) * t174;
t101 = t256 * t153 - t189 * t252;
t308 = t219 * t256;
t76 = pkin(4) * t217 - pkin(10) * t308 + t101;
t309 = t219 * t252;
t87 = -pkin(10) * t309 + t102;
t34 = t251 * t76 + t255 * t87;
t324 = Ifges(5,4) * t252;
t323 = Ifges(5,4) * t256;
t322 = Ifges(5,6) * t252;
t319 = t167 * mrSges(7,3);
t318 = t253 * mrSges(4,1);
t316 = t257 * mrSges(4,2);
t314 = qJ(6) * t218;
t313 = t121 * t358;
t246 = pkin(4) * t299;
t247 = t253 * t321;
t229 = t247 + t246;
t307 = t229 * t174;
t232 = -mrSges(5,1) * t256 + mrSges(5,2) * t252;
t303 = t253 * t232;
t295 = 2 * mrSges(6,3);
t294 = 0.2e1 * mrSges(7,3);
t293 = 0.2e1 * t258;
t290 = t255 * t167 * mrSges(6,3);
t22 = pkin(10) * t311 + pkin(4) * t170 + (-t179 + (pkin(10) * t219 - t153) * t252) * qJD(4) + t280;
t26 = -t265 * pkin(10) + t31;
t8 = -t34 * qJD(5) + t255 * t22 - t251 * t26;
t2 = pkin(5) * t170 - qJ(6) * t48 + qJD(6) * t143 + t8;
t35 = mrSges(7,1) * t170 - mrSges(7,3) * t48;
t289 = m(7) * t2 + t35;
t243 = -pkin(4) * t256 - pkin(3);
t17 = -t49 * mrSges(7,1) + t48 * mrSges(7,2);
t284 = (-mrSges(6,2) - mrSges(7,2)) * t255;
t283 = -t299 / 0.2e1;
t282 = -(2 * Ifges(4,4)) - t322;
t141 = pkin(5) * t168 + t246;
t33 = -t251 * t87 + t255 * t76;
t146 = t255 * t213 - t214 * t251;
t186 = t255 * t235 - t237 * t251;
t275 = mrSges(5,3) * t357;
t274 = m(7) * t55 + t319;
t273 = m(7) * t68 + t319;
t135 = pkin(4) * t309 - t358;
t272 = mrSges(5,1) * t252 + mrSges(5,2) * t256;
t271 = Ifges(5,1) * t256 - t324;
t270 = -Ifges(5,2) * t252 + t323;
t269 = Ifges(5,5) * t252 + Ifges(5,6) * t256;
t193 = pkin(5) * t268 + t243;
t266 = t373 * t170 + t374 * t49 + t375 * t48;
t7 = t251 * t22 + t255 * t26 + t76 * t296 - t87 * t297;
t225 = t270 * qJD(4);
t226 = t271 * qJD(4);
t234 = Ifges(5,1) * t252 + t323;
t263 = -t365 * t167 - t366 * t168 + t367 * t218 + t256 * t225 + t252 * t226 + t234 * t298 - t368 * t268;
t262 = t376 + (mrSges(6,3) + mrSges(7,3)) * (-t168 * t251 + t218 * t297 - t268 * t296) * pkin(4);
t4 = qJ(6) * t49 - qJD(6) * t142 + t7;
t261 = t8 * mrSges(6,1) + t2 * mrSges(7,1) - t7 * mrSges(6,2) - t4 * mrSges(7,2) + t266;
t63 = t265 * pkin(4) + t121;
t154 = -mrSges(5,2) * t217 - mrSges(5,3) * t309;
t155 = mrSges(5,1) * t217 - mrSges(5,3) * t308;
t85 = mrSges(5,1) * t170 + t264 * mrSges(5,3);
t86 = -mrSges(5,2) * t170 - t265 * mrSges(5,3);
t260 = t256 * t86 + m(5) * (-t101 * t298 - t102 * t299 + t364) - t252 * t85 - t155 * t298 - t154 * t299;
t130 = Ifges(5,6) * t217 + t270 * t219;
t131 = Ifges(5,5) * t217 + t271 * t219;
t224 = t272 * qJD(4);
t23 = -pkin(5) * t49 + t63;
t233 = Ifges(5,2) * t256 + t324;
t27 = pkin(5) * t217 + qJ(6) * t143 + t33;
t28 = -qJ(6) * t142 + t34;
t57 = -t264 * Ifges(5,4) - t265 * Ifges(5,2) + Ifges(5,6) * t170;
t58 = -t264 * Ifges(5,1) - t265 * Ifges(5,4) + Ifges(5,5) * t170;
t88 = pkin(5) * t142 + t135;
t259 = ((-t101 * t256 - t102 * t252) * qJD(4) + t364) * mrSges(5,3) - t265 * t233 / 0.2e1 + (-Ifges(5,6) * t299 + t376) * t217 / 0.2e1 + t130 * t283 + (t219 * t283 - t311 / 0.2e1) * t234 + (t167 * t33 - t34 * t168 - t218 * t8 - t268 * t7) * mrSges(6,3) + (t167 * t27 - t28 * t168 - t2 * t218 - t268 * t4) * mrSges(7,3) + t256 * t57 / 0.2e1 + t252 * t58 / 0.2e1 - Ifges(4,5) * t169 - Ifges(4,6) * t170 + t23 * t173 + t63 * t174 + t135 * t104 - t120 * mrSges(4,2) + t88 * t103 - t358 * t224 + t365 * t48 / 0.2e1 + t366 * t49 / 0.2e1 - t367 * t143 / 0.2e1 - t368 * t142 / 0.2e1 - t369 * t167 / 0.2e1 - t370 * t168 / 0.2e1 + t371 * t218 / 0.2e1 + (t232 - mrSges(4,1)) * t121 - t372 * t268 / 0.2e1 + (t375 * t218 - t374 * t268 + t269) * t170 / 0.2e1 + t131 * t298 / 0.2e1 + t226 * t308 / 0.2e1 - t225 * t309 / 0.2e1;
t242 = -pkin(3) - t330;
t241 = pkin(5) + t362;
t231 = t243 - t330;
t211 = t268 * qJ(6);
t192 = t193 - t330;
t150 = t272 * t219;
t138 = -t211 + t188;
t137 = t186 - t314;
t134 = t141 + t247;
t129 = -t211 + t147;
t128 = t146 - t314;
t127 = mrSges(6,1) * t217 + mrSges(6,3) * t143;
t126 = mrSges(7,1) * t217 + mrSges(7,3) * t143;
t125 = -mrSges(6,2) * t217 - mrSges(6,3) * t142;
t124 = -mrSges(7,2) * t217 - mrSges(7,3) * t142;
t90 = mrSges(6,1) * t142 - mrSges(6,2) * t143;
t89 = mrSges(7,1) * t142 - mrSges(7,2) * t143;
t38 = -mrSges(6,2) * t170 + mrSges(6,3) * t49;
t37 = -mrSges(7,2) * t170 + mrSges(7,3) * t49;
t36 = mrSges(6,1) * t170 - mrSges(6,3) * t48;
t18 = -mrSges(6,1) * t49 + mrSges(6,2) * t48;
t1 = [(mrSges(4,1) * t170 - mrSges(4,2) * t169) * t345 + (mrSges(4,3) * t347 - 0.2e1 * Ifges(4,1) * t169 - t252 * t57 + t256 * t58 + (Ifges(5,5) * t256 + t282) * t170 + (-t256 * t130 - t252 * t131 - t217 * t269) * qJD(4)) * t219 + (t101 * t32 + t102 * t31 - t313) * t352 + 0.2e1 * m(4) * (t120 * t189 - t313) + 0.2e1 * t31 * t154 + 0.2e1 * t32 * t155 + 0.2e1 * t135 * t18 + 0.2e1 * t4 * t124 + 0.2e1 * t7 * t125 + 0.2e1 * t2 * t126 + 0.2e1 * t8 * t127 + 0.2e1 * t101 * t85 + 0.2e1 * t102 * t86 + 0.2e1 * t88 * t17 + 0.2e1 * t23 * t89 + 0.2e1 * t63 * t90 + 0.2e1 * t28 * t37 + 0.2e1 * t34 * t38 + 0.2e1 * t27 * t35 + 0.2e1 * t33 * t36 + 0.2e1 * (t169 * t358 - t170 * t189) * mrSges(4,3) - 0.2e1 * t358 * t66 - (t371 + t381) * t143 - (t372 + t382) * t142 + t369 * t48 + t370 * t49 + (-0.2e1 * t120 * mrSges(4,3) - t282 * t169 + ((2 * Ifges(4,2)) + Ifges(5,3) + t373) * t170 + t266 + t359) * t217 - t131 * t311 + t150 * t347 + (t2 * t27 + t23 * t88 + t28 * t4) * t350 + (t135 * t63 + t33 * t8 + t34 * t7) * t351 + t130 * t312 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t258) * t293 + (m(4) * pkin(2) * t345 + 0.2e1 * pkin(2) * (mrSges(4,1) * t217 + mrSges(4,2) * t219) - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t254 + (Ifges(3,1) - Ifges(3,2)) * t293) * t254) * qJD(2); t259 + (Ifges(3,5) * t258 - Ifges(3,6) * t254 + (-mrSges(3,1) * t258 + mrSges(3,2) * t254) * pkin(7)) * qJD(2) + t229 * t90 + t231 * t18 + (m(4) * (t120 * t253 - t121 * t257) + (t169 * t257 - t170 * t253) * mrSges(4,3) + ((-t217 * mrSges(4,3) + t256 * t154 - t252 * t155 + m(5) * (-t101 * t252 + t102 * t256) + m(4) * t189) * t257 + (t219 * mrSges(4,3) + t150 - (m(5) + m(4)) * t358) * t253) * qJD(3)) * pkin(2) + t260 * t240 + t192 * t17 + t146 * t36 + t147 * t38 + t129 * t37 + t134 * t89 + t54 * t124 + t78 * t125 + t55 * t126 + t79 * t127 + t128 * t35 + m(6) * (t135 * t229 + t146 * t8 + t147 * t7 + t231 * t63 + t33 * t79 + t34 * t78) + m(7) * (t128 * t2 + t129 * t4 + t134 * t88 + t192 * t23 + t27 * t55 + t28 * t54) + t377 * t242; -t233 * t299 + 0.2e1 * t242 * t224 + 0.2e1 * t307 + t231 * t348 + t192 * t349 + t134 * t346 + t263 + ((t240 * t357 + t242 * t253) * t352 - 0.2e1 * t316 - 0.2e1 * t318 + 0.2e1 * t303 + 0.2e1 * t275) * t321 + (t128 * t167 - t129 * t168 - t55 * t218 - t268 * t54) * t294 + (t146 * t167 - t147 * t168 - t79 * t218 - t268 * t78) * t295 + (t146 * t79 + t147 * t78 + t229 * t231) * t351 + (t128 * t55 + t129 * t54 + t134 * t192) * t350; t259 + m(6) * (t118 * t34 + t119 * t33 + t135 * t246 + t186 * t8 + t188 * t7 + t243 * t63) + m(7) * (t137 * t2 + t138 * t4 + t141 * t88 + t193 * t23 + t27 * t68 + t28 * t67) + t90 * t246 + t243 * t18 + t260 * pkin(9) + t193 * t17 + t186 * t36 + t188 * t38 + t141 * t89 + t137 * t35 + t138 * t37 + t67 * t124 + t118 * t125 + t68 * t126 + t119 * t127 - t377 * pkin(3); ((-t119 - t79) * t218 - (t118 + t78) * t268 - (t147 + t188) * t168 - (-t146 - t186) * t167) * mrSges(6,3) + ((-t55 - t68) * t218 - (t54 + t67) * t268 - (t129 + t138) * t168 - (-t128 - t137) * t167) * mrSges(7,3) + t307 + m(6) * (t118 * t147 + t119 * t146 + t186 * t79 + t188 * t78 + t229 * t243 + t231 * t246) + m(7) * (t128 * t68 + t129 * t67 + t134 * t193 + t137 * t55 + t138 * t54 + t141 * t192) + (t242 - pkin(3)) * t224 + (t141 + t134) * t173 + t263 + (t243 + t231) * t104 + (t192 + t193) * t103 + (-t233 + t329) * t299 + (m(5) * (-pkin(3) * t253 + pkin(9) * t357) - t316 - t318 + t303 + t275) * t321; (t118 * t188 + t119 * t186 + t243 * t246) * t351 + (t137 * t68 + t138 * t67 + t141 * t193) * t350 + t243 * t348 - 0.2e1 * pkin(3) * t224 + t193 * t349 + t141 * t346 + t263 + (-t233 + 0.2e1 * t329) * t299 + (-t118 * t268 - t119 * t218 + t167 * t186 - t188 * t168) * t295 + (t137 * t167 - t138 * t168 - t68 * t218 - t268 * t67) * t294; -t265 * Ifges(5,6) - Ifges(5,5) * t286 + t289 * t241 + (t255 * t36 + (t37 + t38) * t251 + ((t124 + t125) * t255 + (-t126 - t127) * t251) * qJD(5) + m(7) * (t251 * t4 - t27 * t297 + t28 * t296) + m(6) * (t251 * t7 + t255 * t8 + t34 * t296 - t33 * t297)) * pkin(4) + t261 - t31 * mrSges(5,2) + t32 * mrSges(5,1) + t359; (t232 * t240 - t322) * qJD(4) - t272 * t291 + t274 * t241 + (t290 + m(7) * (-t128 * t297 + t129 * t296 + t251 * t54) + m(6) * (-t146 * t297 + t147 * t296 + t251 * t78 + t255 * t79)) * pkin(4) + t262 + t354; (t232 * pkin(9) - t322) * qJD(4) + t273 * t241 + (t290 + m(7) * (-t137 * t297 + t138 * t296 + t251 * t67) + m(6) * (t118 * t251 + t119 * t255 - t186 * t297 + t188 * t296)) * pkin(4) + t262 + t353; 0.2e1 * (t284 + ((-t241 + t362) * m(7) + t361) * t251) * t360; t289 * pkin(5) + t261; t274 * pkin(5) + t277 + t354; t273 * pkin(5) + t277 + t353; (t284 + (-m(7) * pkin(5) + t361) * t251) * t360; 0; m(7) * t23 + t17; m(7) * t134 + t103; m(7) * t141 + t103; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
