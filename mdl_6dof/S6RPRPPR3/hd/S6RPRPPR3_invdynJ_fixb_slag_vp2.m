% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:07
% EndTime: 2019-03-09 02:44:22
% DurationCPUTime: 11.24s
% Computational Cost: add. (3177->552), mult. (6052->691), div. (0->0), fcn. (3105->10), ass. (0->250)
t366 = Ifges(6,1) + Ifges(5,3);
t347 = Ifges(5,4) + Ifges(4,5);
t365 = -Ifges(4,6) + Ifges(5,6);
t364 = Ifges(6,6) + t347;
t315 = -pkin(4) - pkin(8);
t154 = -pkin(3) + t315;
t168 = -pkin(3) - pkin(4);
t363 = -m(6) * t168 - m(7) * t154 - mrSges(6,2) + mrSges(7,3);
t163 = sin(qJ(3));
t158 = sin(pkin(9));
t119 = pkin(1) * t158 + pkin(7);
t101 = t119 * qJD(1);
t166 = cos(qJ(3));
t263 = qJD(3) * t166;
t356 = qJD(2) * qJD(3) + t119 * qJDD(1);
t27 = qJDD(2) * t166 - t101 * t263 - t163 * t356;
t181 = qJDD(4) - t27;
t256 = qJD(1) * qJD(5);
t257 = qJD(1) * qJD(3);
t98 = qJDD(1) * t163 + t166 * t257;
t171 = -qJ(5) * t98 - t163 * t256 + t181;
t11 = qJDD(3) * t154 + t171;
t162 = sin(qJ(6));
t165 = cos(qJ(6));
t159 = cos(pkin(9));
t281 = pkin(1) * qJDD(1);
t122 = t159 * t281;
t100 = -qJDD(1) * pkin(2) - t122;
t139 = t163 * qJD(4);
t97 = -t166 * qJDD(1) + t163 * t257;
t24 = t97 * pkin(3) - t98 * qJ(4) - qJD(1) * t139 + t100;
t190 = qJDD(5) - t24;
t7 = pkin(5) * t98 + t315 * t97 + t190;
t221 = t163 * pkin(5) + t166 * pkin(8);
t267 = qJD(1) * t166;
t120 = -t159 * pkin(1) - pkin(2);
t105 = t120 * qJD(1);
t259 = t163 * qJD(1);
t53 = -pkin(3) * t267 - qJ(4) * t259 + t105;
t41 = pkin(4) * t267 + qJD(5) - t53;
t25 = qJD(1) * t221 + t41;
t125 = qJ(5) * t259;
t282 = -t166 * qJD(2) + t163 * t101;
t242 = qJD(4) + t282;
t223 = -t125 + t242;
t35 = qJD(3) * t154 + t223;
t8 = -t162 * t35 + t165 * t25;
t1 = qJD(6) * t8 + t11 * t165 + t162 * t7;
t362 = t1 * mrSges(7,2);
t9 = t162 * t25 + t165 * t35;
t2 = -qJD(6) * t9 - t11 * t162 + t165 * t7;
t361 = t2 * mrSges(7,1);
t151 = qJ(1) + pkin(9);
t135 = sin(t151);
t360 = g(2) * t135;
t359 = Ifges(6,5) - t365;
t13 = qJDD(3) * t168 + t171;
t156 = qJD(3) * qJ(4);
t60 = t163 * qJD(2) + t166 * t101;
t45 = -qJ(5) * t267 + t60;
t42 = -t156 - t45;
t358 = -qJD(3) * t42 - t13;
t23 = -qJDD(3) * pkin(3) + t181;
t52 = t156 + t60;
t357 = qJD(3) * t52 - t23;
t214 = t166 * mrSges(5,1) + t163 * mrSges(5,3);
t216 = mrSges(4,1) * t166 - mrSges(4,2) * t163;
t355 = t214 + t216;
t219 = t1 * t165 - t162 * t2;
t261 = qJD(6) * t165;
t262 = qJD(6) * t162;
t354 = t8 * t261 + t9 * t262 - t219;
t145 = t163 * mrSges(6,1);
t339 = -t166 * mrSges(6,2) + t145;
t353 = -t166 * mrSges(7,3) - t339;
t264 = qJD(3) * t165;
t92 = t162 * t267 - t264;
t33 = qJD(6) * t92 - qJDD(3) * t162 + t165 * t97;
t352 = -t33 / 0.2e1;
t266 = qJD(3) * t162;
t93 = t165 * t267 + t266;
t34 = qJD(6) * t93 - qJDD(3) * t165 - t162 * t97;
t351 = -t34 / 0.2e1;
t86 = qJDD(6) + t98;
t350 = -t86 / 0.2e1;
t349 = -m(7) - m(6);
t218 = t9 * t162 + t8 * t165;
t348 = m(7) * (-qJD(6) * t218 + t219);
t10 = -mrSges(7,1) * t34 + mrSges(7,2) * t33;
t62 = -qJDD(3) * mrSges(6,1) - mrSges(6,3) * t97;
t345 = t10 - t62;
t67 = -mrSges(5,2) * t97 + qJDD(3) * mrSges(5,3);
t344 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t97 + t67;
t211 = mrSges(7,1) * t162 + mrSges(7,2) * t165;
t38 = qJD(3) * pkin(5) - t42;
t343 = t211 * t38;
t243 = mrSges(6,3) * t267;
t283 = qJD(3) * mrSges(6,1) - mrSges(7,1) * t92 - mrSges(7,2) * t93 - t243;
t304 = Ifges(4,4) * t163;
t206 = t166 * Ifges(4,2) + t304;
t115 = qJD(6) + t259;
t307 = t93 * Ifges(7,4);
t29 = t92 * Ifges(7,2) + t115 * Ifges(7,6) - t307;
t341 = Ifges(4,6) * qJD(3) + qJD(1) * t206 + t162 * t29;
t246 = mrSges(4,3) * t259;
t248 = mrSges(5,2) * t259;
t340 = t246 + t248 + (-mrSges(4,1) - mrSges(5,1)) * qJD(3);
t338 = t163 * t365 + t166 * t347;
t18 = mrSges(7,1) * t86 - mrSges(7,3) * t33;
t19 = -mrSges(7,2) * t86 + mrSges(7,3) * t34;
t197 = t162 * t18 - t165 * t19;
t265 = qJD(3) * t163;
t26 = t163 * qJDD(2) - t101 * t265 + t166 * t356;
t336 = -t163 * t27 + t166 * t26;
t22 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t26;
t335 = t163 * t23 + t166 * t22;
t136 = cos(t151);
t334 = g(1) * t136 + t360;
t258 = m(5) - t349;
t132 = Ifges(4,4) * t267;
t297 = Ifges(5,5) * t166;
t210 = t163 * Ifges(5,1) - t297;
t333 = Ifges(4,1) * t259 - t93 * Ifges(7,5) + t92 * Ifges(7,6) + t115 * Ifges(7,3) + qJD(1) * t210 + qJD(3) * t347 + t132;
t131 = Ifges(5,5) * t259;
t302 = Ifges(6,4) * t163;
t209 = -t166 * Ifges(6,1) - t302;
t80 = Ifges(7,4) * t92;
t30 = -Ifges(7,1) * t93 + Ifges(7,5) * t115 + t80;
t332 = -Ifges(5,3) * t267 + qJD(1) * t209 + t165 * t30 + t131 + (-Ifges(6,5) + Ifges(5,6)) * qJD(3);
t331 = mrSges(3,1) + t355;
t212 = mrSges(7,1) * t165 - mrSges(7,2) * t162;
t288 = t166 * mrSges(5,3);
t290 = t166 * mrSges(6,1);
t330 = -(m(7) * pkin(5) + t212) * t166 - t290 - t288 + (m(5) * pkin(3) + mrSges(5,1) + t363) * t163;
t301 = Ifges(6,4) * t166;
t329 = t163 * t302 + (t297 - t301 + (-Ifges(6,2) + t366) * t163) * t166;
t247 = mrSges(5,2) * t267;
t108 = qJD(3) * mrSges(5,3) + t247;
t328 = -m(5) * t52 - t108;
t327 = mrSges(3,2) - mrSges(4,3) - mrSges(5,2) + mrSges(6,3);
t326 = qJ(4) * t258;
t46 = -mrSges(7,2) * t115 + mrSges(7,3) * t92;
t47 = mrSges(7,1) * t115 + mrSges(7,3) * t93;
t195 = t162 * t46 + t165 * t47;
t65 = -qJDD(3) * mrSges(5,1) + t98 * mrSges(5,2);
t66 = qJDD(3) * mrSges(6,2) - t98 * mrSges(6,3);
t325 = t195 * qJD(6) + t197 - t65 - t66;
t324 = m(6) * t42 - m(7) * t38 - t283;
t321 = Ifges(7,1) * t352 + Ifges(7,4) * t351 + Ifges(7,5) * t350;
t319 = -t93 / 0.2e1;
t318 = t93 / 0.2e1;
t164 = sin(qJ(1));
t312 = pkin(1) * t164;
t149 = t166 * pkin(3);
t167 = cos(qJ(1));
t150 = t167 * pkin(1);
t306 = -qJD(1) / 0.2e1;
t305 = -qJD(6) / 0.2e1;
t161 = qJ(4) + pkin(5);
t303 = Ifges(4,4) * t166;
t300 = Ifges(7,4) * t162;
t299 = Ifges(7,4) * t165;
t298 = Ifges(5,5) * t163;
t280 = qJD(3) * t38;
t276 = t119 * t166;
t275 = t136 * t166;
t274 = t162 * t163;
t273 = t162 * t166;
t142 = t163 * qJ(4);
t272 = t163 * t165;
t271 = t165 * t166;
t270 = -qJ(5) + t119;
t269 = qJ(4) * t263 + t139;
t268 = t149 + t142;
t260 = qJD(6) * t166;
t252 = Ifges(7,5) * t33 + Ifges(7,6) * t34 + Ifges(7,3) * t86;
t249 = -t108 - t283;
t245 = mrSges(4,3) * t267;
t244 = mrSges(6,3) * t259;
t241 = t136 * pkin(2) + t135 * pkin(7) + t150;
t148 = t166 * pkin(4);
t240 = t148 + t268;
t239 = t168 * qJD(3);
t238 = t119 * t265;
t236 = t98 * mrSges(6,1) + t97 * mrSges(6,2);
t228 = t136 * pkin(7) - t312;
t227 = -pkin(2) - t142;
t226 = -t257 / 0.2e1;
t225 = t257 / 0.2e1;
t76 = -t268 + t120;
t220 = g(1) * t135 - g(2) * t136;
t215 = mrSges(4,1) * t163 + mrSges(4,2) * t166;
t208 = Ifges(7,1) * t165 - t300;
t207 = Ifges(7,1) * t162 + t299;
t204 = -t163 * Ifges(6,2) - t301;
t203 = -Ifges(7,2) * t162 + t299;
t202 = Ifges(7,2) * t165 + t300;
t200 = Ifges(6,5) * t163 - Ifges(6,6) * t166;
t199 = Ifges(7,5) * t165 - Ifges(7,6) * t162;
t198 = Ifges(7,5) * t162 + Ifges(7,6) * t165;
t193 = t240 + t221;
t43 = t193 - t120;
t78 = t270 * t163;
t20 = -t162 * t78 + t165 * t43;
t21 = t162 * t43 + t165 * t78;
t196 = -t162 * t47 + t165 * t46;
t194 = pkin(3) * t275 + t136 * t142 + t241;
t102 = qJD(3) * mrSges(6,2) - t244;
t191 = t102 + t196;
t189 = t41 * (t163 * mrSges(6,2) + t290);
t188 = t53 * (mrSges(5,1) * t163 - t288);
t187 = t105 * t215;
t186 = t163 * (Ifges(4,1) * t166 - t304);
t182 = pkin(5) * t166 + t154 * t163;
t179 = t162 * t260 + t163 * t264;
t178 = t162 * t265 - t165 * t260;
t175 = Ifges(7,5) * t166 + t163 * t208;
t174 = Ifges(7,6) * t166 + t163 * t203;
t173 = Ifges(7,3) * t166 + t163 * t199;
t14 = -qJ(5) * t97 + t166 * t256 - t22;
t128 = qJ(4) * t267;
t107 = -qJD(3) * mrSges(4,2) + t245;
t96 = t339 * qJD(1);
t95 = pkin(3) * t259 - t128;
t94 = t214 * qJD(1);
t77 = t211 * t166;
t70 = -Ifges(6,6) * qJD(3) + qJD(1) * t204;
t68 = pkin(3) * t265 - t269;
t64 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t98;
t61 = t168 * t259 + t128;
t58 = t148 - t76;
t57 = -t135 * t162 + t136 * t272;
t56 = -t135 * t165 - t136 * t274;
t55 = -t135 * t272 - t136 * t162;
t54 = t135 * t274 - t136 * t165;
t51 = -qJD(5) * t163 + t263 * t270;
t49 = t163 * t239 + t269;
t48 = -qJD(3) * pkin(3) + t242;
t44 = -t125 + t282;
t40 = qJD(1) * t182 + t128;
t37 = t239 + t223;
t36 = qJD(3) * t182 + t269;
t17 = -pkin(4) * t97 + t190;
t16 = t162 * t40 + t165 * t45;
t15 = -t162 * t45 + t165 * t40;
t12 = qJDD(3) * pkin(5) - t14;
t5 = t33 * Ifges(7,4) + t34 * Ifges(7,2) + t86 * Ifges(7,6);
t4 = -qJD(6) * t21 - t162 * t51 + t165 * t36;
t3 = qJD(6) * t20 + t162 * t36 + t165 * t51;
t6 = [m(7) * (t1 * t21 + t2 * t20 + t3 * t9 + t4 * t8) + (-t108 - t107) * t238 + (-t64 + t65) * t119 * t163 + (t252 + (Ifges(4,1) + Ifges(5,1)) * t98 + (-Ifges(4,4) + Ifges(5,5)) * t97) * t163 / 0.2e1 + (-m(3) * t150 - m(4) * t241 - m(5) * t194 - mrSges(2,1) * t167 - t57 * mrSges(7,1) + mrSges(2,2) * t164 - t56 * mrSges(7,2) + t349 * (pkin(4) * t275 - qJ(5) * t135 + t194) + t327 * t135 + (-m(7) * t221 - t331 + t353) * t136) * g(2) - t163 * t362 + ((t359 / 0.2e1 + Ifges(4,6) / 0.2e1 + Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1) * t166 + t364 * t163) * qJDD(3) + (m(3) * t312 + mrSges(2,1) * t164 - t55 * mrSges(7,1) + mrSges(2,2) * t167 - t54 * mrSges(7,2) + (-m(5) - m(4)) * t228 + t349 * (-qJ(5) * t136 + t228) + t327 * t136 + (-m(5) * (t227 - t149) - m(7) * (-t161 * t163 - pkin(2)) - m(6) * t227 + t145 + m(4) * pkin(2) + t363 * t166 + t331) * t135) * g(1) + (-t200 / 0.2e1 + t338 / 0.2e1) * qJD(3) ^ 2 + (-t37 * mrSges(6,3) - t70 / 0.2e1 - t9 * mrSges(7,2) + t8 * mrSges(7,1) + t340 * t119 + t333 / 0.2e1 + t48 * mrSges(5,2) + t282 * mrSges(4,3)) * t263 + m(4) * (t100 * t120 + ((-t60 * t163 + t166 * t282) * qJD(3) + t336) * t119) - t100 * t216 + t17 * t339 + (t1 * t273 - t178 * t9 - t179 * t8 + t2 * t271) * mrSges(7,3) - t163 * (Ifges(6,4) * t97 - Ifges(6,2) * t98) / 0.2e1 - ((-Ifges(6,4) + Ifges(5,5)) * t98 + t366 * t97) * t166 / 0.2e1 + t166 * (Ifges(4,4) * t98 - Ifges(4,2) * t97) / 0.2e1 + (-t42 * mrSges(6,3) - t324 * qJ(5) + t332 / 0.2e1 - t52 * mrSges(5,2) - t60 * mrSges(4,3) - t341 / 0.2e1) * t265 + (m(3) * (t158 ^ 2 + t159 ^ 2) * pkin(1) ^ 2 + Ifges(3,3) + Ifges(2,3)) * qJDD(1) - 0.2e1 * t158 * mrSges(3,2) * t281 + m(6) * (t13 * t78 + t17 * t58 + t37 * t51 + t41 * t49) + t5 * t273 / 0.2e1 + t92 * (qJD(3) * t174 + t202 * t260) / 0.2e1 + t115 * (qJD(3) * t173 + t198 * t260) / 0.2e1 + 0.2e1 * mrSges(3,1) * t122 + t58 * t236 + (-t13 * t163 + t14 * t166) * mrSges(6,3) + t335 * mrSges(5,2) + t336 * mrSges(4,3) + t120 * (mrSges(4,1) * t97 + mrSges(4,2) * t98) + t76 * (mrSges(5,1) * t97 - mrSges(5,3) * t98) + t51 * t102 - t68 * t94 + t49 * t96 - t12 * t77 + t78 * t66 + t3 * t46 + t4 * t47 + t20 * t18 + t21 * t19 + (t162 * t30 + t165 * t29) * t260 / 0.2e1 + (t163 * Ifges(4,1) + t210 + t303) * t98 / 0.2e1 + (-t166 * Ifges(5,3) + t209 + t298) * t97 / 0.2e1 + t324 * (qJD(5) * t166 + t238) + (qJD(3) * t175 + t207 * t260) * t319 + t271 * t321 + (t163 * (Ifges(5,1) * t166 + t298) + t166 * (-Ifges(4,2) * t163 + t303) + t186) * t225 - t24 * t214 + t34 * (Ifges(7,6) * t163 - t166 * t203) / 0.2e1 - t98 * t204 / 0.2e1 - t97 * t206 / 0.2e1 + t33 * (Ifges(7,5) * t163 - t166 * t208) / 0.2e1 + t86 * (Ifges(7,3) * t163 - t166 * t199) / 0.2e1 + t163 * t361 + qJD(3) * t187 + qJD(3) * t188 + qJD(3) * t189 + t38 * (mrSges(7,1) * t178 + mrSges(7,2) * t179) + t329 * t226 + m(5) * (t24 * t76 + t53 * t68 + ((-t52 * t163 + t48 * t166) * qJD(3) + t335) * t119) + t344 * t276 + (-m(6) * t14 + m(7) * t12 + t345) * (-t166 * qJ(5) + t276); m(3) * qJDD(2) + (-m(3) - m(4) - t258) * g(3) + ((t191 + t340) * qJD(3) + m(7) * (t264 * t9 - t266 * t8 + t12) + m(5) * (qJD(3) * t48 + t22) + m(4) * (qJD(3) * t282 + t26) + m(6) * (qJD(3) * t37 - t14) + t344 + t345) * t163 + (t64 + (t107 - t249) * qJD(3) + m(7) * (t280 + t354) + m(5) * t357 + m(4) * (qJD(3) * t60 + t27) + m(6) * t358 + t325) * t166; t354 * mrSges(7,3) + (-m(5) * t268 - m(6) * t240 - m(7) * t193 - t163 * t212 + t353 - t355) * g(3) + (-pkin(3) * t23 + qJ(4) * t22 - t53 * t95) * m(5) + t364 * t98 + (t12 * t161 - t15 * t8 - t16 * t9 + t38 * t44) * m(7) + (-qJ(4) * t14 + t13 * t168 - t37 * t45 - t41 * t61 - t42 * t44) * m(6) + (t174 * t306 + t203 * t305) * t92 - (-t107 + t245 + t328) * t282 + (t136 * t330 - t275 * t326) * g(1) - t359 * t97 + (Ifges(6,3) + Ifges(4,3) + Ifges(5,2)) * qJDD(3) + t70 * t267 / 0.2e1 - t30 * t261 / 0.2e1 + t29 * t262 / 0.2e1 - t48 * t247 + t168 * t66 - t165 * t5 / 0.2e1 + t161 * t10 - t45 * t102 + t95 * t94 - t61 * t96 - pkin(3) * t65 - t16 * t46 - t15 * t47 + t22 * mrSges(5,3) - t23 * mrSges(5,1) - t26 * mrSges(4,2) + t27 * mrSges(4,1) + t198 * t350 + t202 * t351 + t207 * t352 + (t173 * t306 + t199 * t305) * t115 + (t67 - t62) * qJ(4) + t13 * mrSges(6,2) - t14 * mrSges(6,1) + t37 * t243 + t42 * t244 + t52 * t248 + t162 * t321 + (-t9 * (-mrSges(7,2) * t166 - mrSges(7,3) * t274) - t8 * (mrSges(7,1) * t166 - mrSges(7,3) * t272) + t175 * t318 - t187 - t188 - t189 + (-t186 / 0.2e1 + t329 / 0.2e1) * qJD(1)) * qJD(1) + t12 * t212 + t200 * t225 + (-t166 * t326 + t330) * t360 + (-t324 - t328) * qJD(4) - (Ifges(5,1) * t267 + t131 + t332) * t259 / 0.2e1 - (-Ifges(4,2) * t259 + t132 + t333) * t267 / 0.2e1 + t334 * t215 + t338 * t226 + (-m(5) * t48 + t246 - t340) * t60 + t341 * t259 / 0.2e1 + t283 * t44 - t259 * t343 + (t208 * t318 - t343) * qJD(6) + (-t261 * t47 - t262 * t46 - t197 + t348) * t154; t249 * qJD(3) + t258 * t166 * g(3) + ((-t195 - t94 - t96) * qJD(1) - t334 * t258) * t163 - m(7) * (t218 * t259 + t280) + t348 - t325 + (-t259 * t41 - t358) * m(6) + (t259 * t53 - t357) * m(5); t162 * t19 + t165 * t18 + t196 * qJD(6) + (t1 * t162 + t165 * t2 + (-t162 * t8 + t165 * t9) * qJD(6) + t220) * m(7) + (t17 + t220) * m(6) + (t283 * t166 + t191 * t163 - m(6) * (-t163 * t37 + t166 * t42) - m(7) * (-t166 * t38 - t272 * t9 + t274 * t8)) * qJD(1) + t236; -t362 + t361 - t38 * (-mrSges(7,1) * t93 + mrSges(7,2) * t92) + (Ifges(7,1) * t92 + t307) * t318 + t29 * t319 - t115 * (Ifges(7,5) * t92 + Ifges(7,6) * t93) / 0.2e1 - t8 * t46 + t9 * t47 - g(1) * (mrSges(7,1) * t56 - mrSges(7,2) * t57) - g(2) * (-mrSges(7,1) * t54 + mrSges(7,2) * t55) - g(3) * t77 + (t8 * t92 - t9 * t93) * mrSges(7,3) + t252 - (Ifges(7,2) * t93 + t30 + t80) * t92 / 0.2e1;];
tau  = t6;
