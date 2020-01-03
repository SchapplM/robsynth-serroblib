% Calculate vector of inverse dynamics joint torques for
% S5RRRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:00
% EndTime: 2019-12-31 21:11:11
% DurationCPUTime: 5.50s
% Computational Cost: add. (3493->416), mult. (5202->543), div. (0->0), fcn. (2801->10), ass. (0->207)
t376 = mrSges(4,1) + mrSges(5,1);
t371 = Ifges(4,5) + Ifges(5,4);
t370 = Ifges(5,6) - Ifges(4,6);
t200 = qJD(1) + qJD(2);
t209 = sin(qJ(2));
t318 = pkin(1) * qJD(1);
t286 = t209 * t318;
t148 = pkin(7) * t200 + t286;
t208 = sin(qJ(3));
t129 = t208 * t148;
t375 = qJD(4) + t129;
t199 = -qJD(3) + qJD(5);
t374 = t199 / 0.2e1;
t298 = t200 * t208;
t373 = -pkin(8) * t298 + t375;
t281 = mrSges(5,2) * t298;
t354 = -mrSges(4,3) * t298 + t376 * qJD(3) - t281;
t212 = cos(qJ(3));
t297 = t200 * t212;
t147 = mrSges(5,2) * t297 + qJD(3) * mrSges(5,3);
t353 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t297 + t147;
t372 = -mrSges(5,2) + mrSges(6,3) - mrSges(4,3) + mrSges(3,2);
t170 = Ifges(4,4) * t297;
t320 = Ifges(5,5) * t212;
t245 = Ifges(5,1) * t208 - t320;
t369 = Ifges(4,1) * t298 + t371 * qJD(3) + t200 * t245 + t170;
t368 = t370 * t208 + t371 * t212;
t213 = cos(qJ(2));
t317 = pkin(1) * qJD(2);
t279 = qJD(1) * t317;
t306 = pkin(1) * qJDD(1);
t137 = t209 * t306 + t213 * t279;
t198 = qJDD(1) + qJDD(2);
t113 = pkin(7) * t198 + t137;
t289 = qJD(3) * t208;
t47 = t212 * t113 - t148 * t289;
t288 = qJD(3) * t212;
t48 = -t208 * t113 - t148 * t288;
t367 = -t208 * t48 + t212 * t47;
t38 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t47;
t270 = qJDD(4) - t48;
t43 = -qJDD(3) * pkin(3) + t270;
t366 = t208 * t43 + t212 * t38;
t285 = t213 * t318;
t149 = -pkin(2) * t200 - t285;
t248 = mrSges(4,1) * t208 + mrSges(4,2) * t212;
t105 = -qJD(3) * pkin(3) + t375;
t304 = t105 * t212;
t310 = t212 * mrSges(5,3);
t195 = t212 * pkin(3);
t193 = t208 * qJ(4);
t266 = pkin(2) + t193;
t232 = -t266 - t195;
t66 = t200 * t232 - t285;
t365 = t149 * t248 + t66 * (t208 * mrSges(5,1) - t310) + mrSges(5,2) * t304;
t123 = -t212 * t198 + t200 * t289;
t124 = t198 * t208 + t200 * t288;
t130 = t212 * t148;
t203 = qJD(3) * qJ(4);
t115 = t130 + t203;
t237 = -t115 * t208 + t304;
t93 = -qJDD(3) * mrSges(5,1) + t124 * mrSges(5,2);
t94 = -mrSges(5,2) * t123 + qJDD(3) * mrSges(5,3);
t364 = t208 * (-qJDD(3) * mrSges(4,1) + mrSges(4,3) * t124 + t93) + t212 * (-qJDD(3) * mrSges(4,2) - mrSges(4,3) * t123 + t94) + m(4) * t367 + m(5) * (qJD(3) * t237 + t366) - t353 * t289 - t354 * t288;
t207 = sin(qJ(5));
t211 = cos(qJ(5));
t352 = -t207 * t212 + t208 * t211;
t111 = t352 * t200;
t100 = Ifges(6,4) * t111;
t235 = t207 * t208 + t211 * t212;
t110 = t235 * t200;
t197 = -qJDD(3) + qJDD(5);
t339 = pkin(3) + pkin(4);
t23 = -pkin(8) * t124 - qJDD(3) * t339 + t270;
t55 = -qJD(3) * t339 + t373;
t84 = -pkin(8) * t297 + t130;
t67 = t203 + t84;
t24 = -t207 * t67 + t211 * t55;
t26 = pkin(8) * t123 + t38;
t2 = qJD(5) * t24 + t207 * t23 + t211 * t26;
t25 = t207 * t55 + t211 * t67;
t225 = t235 * qJD(5);
t28 = t123 * t207 + t124 * t211 - t200 * t225;
t226 = t352 * qJD(5);
t29 = t123 * t211 - t124 * t207 - t200 * t226;
t3 = -qJD(5) * t25 - t207 * t26 + t211 * t23;
t335 = t111 / 0.2e1;
t35 = -Ifges(6,2) * t110 + Ifges(6,6) * t199 + t100;
t99 = Ifges(6,4) * t110;
t36 = Ifges(6,1) * t111 + Ifges(6,5) * t199 - t99;
t346 = t212 * t339 + t266;
t54 = t200 * t346 + t285;
t363 = -(-t110 * t24 + t111 * t25) * mrSges(6,3) - t111 * (Ifges(6,1) * t110 + t100) / 0.2e1 - t3 * mrSges(6,1) + t2 * mrSges(6,2) - Ifges(6,5) * t28 - Ifges(6,6) * t29 - Ifges(6,3) * t197 - t335 * t35 + (-Ifges(6,5) * t110 - Ifges(6,6) * t111) * t374 + t54 * (mrSges(6,1) * t111 - mrSges(6,2) * t110) + (Ifges(6,2) * t111 - t36 + t99) * t110 / 0.2e1;
t362 = m(5) + m(6);
t187 = pkin(8) * t289;
t142 = -pkin(7) * t289 + t187;
t338 = pkin(7) - pkin(8);
t159 = t338 * t212;
t143 = qJD(3) * t159;
t158 = t338 * t208;
t73 = t158 * t207 + t159 * t211;
t361 = -qJD(5) * t73 - t142 * t207 + t143 * t211 - t285 * t352;
t72 = t158 * t211 - t159 * t207;
t360 = qJD(5) * t72 + t142 * t211 + t143 * t207 - t235 * t285;
t150 = -t207 * qJ(4) - t211 * t339;
t359 = qJD(5) * t150 - t207 * t84 + t373 * t211;
t151 = t211 * qJ(4) - t207 * t339;
t358 = -qJD(5) * t151 - t373 * t207 - t211 * t84;
t247 = t212 * mrSges(5,1) + t208 * mrSges(5,3);
t120 = t247 * t200;
t44 = mrSges(6,1) * t110 + mrSges(6,2) * t111;
t357 = -t120 - t44;
t263 = mrSges(6,1) * t235 + mrSges(6,2) * t352;
t206 = qJ(1) + qJ(2);
t191 = sin(t206);
t192 = cos(t206);
t350 = g(1) * t192 + g(2) * t191;
t136 = -t209 * t279 + t213 * t306;
t253 = pkin(2) * t198 + t136;
t349 = m(4) * t253 - mrSges(4,1) * t123 - mrSges(4,2) * t124;
t301 = t191 * t212;
t313 = t208 * mrSges(4,2);
t95 = t352 * t191;
t96 = t235 * t191;
t345 = mrSges(4,1) * t301 + t96 * mrSges(6,1) + t95 * mrSges(6,2) + t372 * t192 + (-m(5) * t232 + m(6) * t346 + mrSges(3,1) + t247 - t313) * t191;
t299 = t192 * t212;
t97 = t352 * t192;
t98 = t235 * t192;
t344 = -t98 * mrSges(6,1) - t97 * mrSges(6,2) - t376 * t299 + t372 * t191 + (-mrSges(3,1) + (mrSges(4,2) - mrSges(5,3)) * t208) * t192;
t330 = pkin(1) * t213;
t332 = pkin(1) * t209;
t343 = m(4) * pkin(1) * (t149 * t209 + (t208 ^ 2 + t212 ^ 2) * t213 * t148) + (-mrSges(3,1) * t332 - mrSges(3,2) * t330) * t200;
t333 = t208 / 0.2e1;
t210 = sin(qJ(1));
t331 = pkin(1) * t210;
t214 = cos(qJ(1));
t196 = t214 * pkin(1);
t183 = pkin(7) + t332;
t325 = -pkin(8) + t183;
t323 = Ifges(4,4) * t208;
t322 = Ifges(4,4) * t212;
t321 = Ifges(5,5) * t208;
t305 = qJ(4) * t212;
t292 = t192 * pkin(2) + t191 * pkin(7);
t291 = t195 + t193;
t290 = qJD(3) * t200;
t190 = t208 * qJD(4);
t284 = t209 * t317;
t283 = t213 * t317;
t280 = t339 * t208;
t119 = pkin(3) * t289 - qJ(4) * t288 - t190;
t194 = t212 * pkin(4);
t277 = t194 + t291;
t184 = -pkin(2) - t330;
t135 = t325 * t212;
t269 = -t290 / 0.2e1;
t180 = t192 * pkin(7);
t265 = -pkin(2) * t191 + t180;
t264 = -pkin(8) * t192 + t180;
t257 = t208 * t283;
t256 = t212 * t283;
t255 = pkin(3) * t299 + t192 * t193 + t292;
t252 = t95 * mrSges(6,1) - t96 * mrSges(6,2);
t251 = t97 * mrSges(6,1) - t98 * mrSges(6,2);
t154 = -t212 * mrSges(4,1) + t313;
t244 = t212 * Ifges(4,2) + t323;
t64 = -mrSges(6,2) * t199 - t110 * mrSges(6,3);
t65 = mrSges(6,1) * t199 - mrSges(6,3) * t111;
t239 = -t207 * t65 + t211 * t64;
t131 = t184 - t291;
t134 = t325 * t208;
t51 = t134 * t211 - t135 * t207;
t52 = t134 * t207 + t135 * t211;
t229 = t208 * (Ifges(4,1) * t212 - t323);
t228 = t212 * (Ifges(5,3) * t208 + t320);
t85 = -pkin(4) * t289 - t119;
t227 = (t105 * t208 + t115 * t212) * t213;
t224 = pkin(4) * t299 - pkin(8) * t191 + t255;
t223 = t310 + (-m(5) * pkin(3) - mrSges(5,1)) * t208;
t221 = qJ(4) * t124 + t190 * t200 + t253;
t169 = Ifges(5,5) * t298;
t106 = Ifges(5,6) * qJD(3) - Ifges(5,3) * t297 + t169;
t107 = Ifges(4,6) * qJD(3) + t200 * t244;
t18 = -t123 * t339 + t221;
t27 = pkin(3) * t123 - t221;
t56 = qJD(3) * t352 - t226;
t57 = qJD(3) * t235 - t225;
t216 = (t371 * t208 - t370 * t212) * qJDD(3) / 0.2e1 + ((Ifges(4,1) + Ifges(5,1)) * t124 + t371 * qJDD(3)) * t333 + (-t115 * t289 + t366) * mrSges(5,2) + t367 * mrSges(4,3) + (t200 * (-Ifges(4,2) * t208 + t322) + t369) * t288 / 0.2e1 - (t3 * mrSges(6,3) - Ifges(6,1) * t28 - Ifges(6,4) * t29 - Ifges(6,5) * t197) * t352 + (t200 * (Ifges(5,1) * t212 + t321) + t106) * t289 / 0.2e1 + (Ifges(4,1) * t208 + t245 + t322) * t124 / 0.2e1 - t253 * t154 + (-t2 * mrSges(6,3) - Ifges(6,4) * t28 - Ifges(6,2) * t29 - Ifges(6,6) * t197) * t235 + (-t244 / 0.2e1 + t321 / 0.2e1 + (Ifges(5,5) - Ifges(4,4)) * t333 + (-Ifges(4,2) / 0.2e1 - Ifges(5,3)) * t212) * t123 + (t365 + t368 * qJD(3) / 0.2e1) * qJD(3) + (Ifges(6,5) * t57 + Ifges(6,6) * t56) * t374 + t212 * (Ifges(4,4) * t124 + Ifges(4,6) * qJDD(3)) / 0.2e1 + Ifges(3,3) * t198 + t136 * mrSges(3,1) - t137 * mrSges(3,2) - t110 * (Ifges(6,4) * t57 + Ifges(6,2) * t56) / 0.2e1 + t228 * t269 + t18 * t263 + t56 * t35 / 0.2e1 + t57 * t36 / 0.2e1 + t54 * (-mrSges(6,1) * t56 + mrSges(6,2) * t57) + (-t24 * t57 + t25 * t56) * mrSges(6,3) - t27 * t247 - t212 * (Ifges(5,5) * t124 + Ifges(5,6) * qJDD(3)) / 0.2e1 - t107 * t289 / 0.2e1 + t229 * t290 / 0.2e1 + (Ifges(6,1) * t57 + Ifges(6,4) * t56) * t335;
t162 = qJ(4) * t299;
t160 = qJ(4) * t301;
t153 = -pkin(2) - t291;
t132 = pkin(2) + t277;
t122 = (pkin(3) * t208 - t305) * t200;
t121 = t154 * t200;
t114 = -t131 + t194;
t86 = t119 + t284;
t78 = (-t280 + t305) * t200;
t76 = qJD(3) * t135 + t257;
t75 = -t183 * t289 + t187 + t256;
t59 = t85 - t284;
t49 = mrSges(5,1) * t123 - mrSges(5,3) * t124;
t22 = -mrSges(6,2) * t197 + mrSges(6,3) * t29;
t21 = mrSges(6,1) * t197 - mrSges(6,3) * t28;
t11 = -qJD(5) * t52 - t207 * t75 + t211 * t76;
t10 = qJD(5) * t51 + t207 * t76 + t211 * t75;
t6 = -mrSges(6,1) * t29 + mrSges(6,2) * t28;
t1 = [m(3) * (t136 * t213 + t137 * t209) * pkin(1) + m(5) * (t131 * t27 + t227 * t317 + t66 * t86) + t131 * t49 + t114 * t6 - t86 * t120 + t59 * t44 + t10 * t64 + t11 * t65 + t51 * t21 + t52 * t22 + t121 * t284 + t216 + m(6) * (t10 * t25 + t11 * t24 + t114 * t18 + t2 * t52 + t3 * t51 + t54 * t59) + Ifges(2,3) * qJDD(1) - t354 * t257 + t353 * t256 + (mrSges(3,1) * t330 - mrSges(3,2) * t332) * t198 - t349 * t184 + t343 * qJD(2) + (-mrSges(2,1) * t214 + mrSges(2,2) * t210 - m(5) * (t196 + t255) - m(6) * (t196 + t224) - m(4) * (t196 + t292) - m(3) * t196 + t344) * g(2) + (mrSges(2,1) * t210 + mrSges(2,2) * t214 + m(3) * t331 - m(5) * (t180 - t331) - m(6) * (t264 - t331) - m(4) * (t265 - t331) + t345) * g(1) + t364 * t183; -t119 * t120 + t132 * t6 + t153 * t49 + t72 * t21 + t73 * t22 + t85 * t44 + t216 + t361 * t65 + t360 * t64 + t349 * pkin(2) + (t132 * t18 + t2 * t73 + t24 * t361 + t25 * t360 + t3 * t72 + t54 * t85) * m(6) + (-(t209 * t66 + t227) * t318 + t119 * t66 + t153 * t27) * m(5) + (m(6) * t54 - t121 - t357) * t286 - t343 * qJD(1) + (-m(4) * t292 - m(5) * t255 - m(6) * t224 + t344) * g(2) + (-m(4) * t265 - m(5) * t180 - m(6) * t264 + t345) * g(1) + (t208 * t354 - t212 * t353) * t285 + t364 * pkin(7); t370 * t123 + t371 * t124 + t368 * t269 - (-Ifges(4,2) * t298 + t170 + t369) * t297 / 0.2e1 + ((-t229 / 0.2e1 + t228 / 0.2e1) * t200 - t365) * t200 + t358 * t65 + (t150 * t3 + t151 * t2 - t54 * t78 - g(1) * (-t192 * t280 + t162) - g(2) * (-t191 * t280 + t160) + t359 * t25 + t358 * t24) * m(6) + t359 * t64 + (Ifges(4,3) + Ifges(5,2)) * qJDD(3) + (-m(5) * t291 - m(6) * t277 + t154 - t247 - t263) * g(3) + (-pkin(3) * t43 - g(1) * t162 - g(2) * t160 + qJ(4) * t38 + qJD(4) * t115 - t122 * t66 - t148 * t237) * m(5) + qJD(4) * t147 + t150 * t21 + t151 * t22 + t122 * t120 - pkin(3) * t93 + qJ(4) * t94 - t78 * t44 - (Ifges(5,1) * t297 + t106 + t169) * t298 / 0.2e1 - t47 * mrSges(4,2) + t48 * mrSges(4,1) + t38 * mrSges(5,3) - t43 * mrSges(5,1) + t115 * t281 + (-t191 * t223 + t252) * g(2) + (-t192 * t223 + t251) * g(1) + t363 + t353 * t129 + t354 * t130 + t350 * t248 + t107 * t298 / 0.2e1; t207 * t22 + t211 * t21 + t239 * qJD(5) + t362 * t212 * g(3) + (-t147 - t239) * qJD(3) + (t357 * t200 - t350 * t362) * t208 + t93 + (t2 * t207 + t211 * t3 - t298 * t54 + t199 * (-t207 * t24 + t211 * t25)) * m(6) + (-qJD(3) * t115 + t298 * t66 + t43) * m(5); -g(1) * t251 - g(2) * t252 + g(3) * t263 - t24 * t64 + t25 * t65 - t363;];
tau = t1;
