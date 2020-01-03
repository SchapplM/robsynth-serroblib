% Calculate vector of inverse dynamics joint torques for
% S5RRRRP2
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:05
% EndTime: 2020-01-03 12:11:15
% DurationCPUTime: 4.42s
% Computational Cost: add. (4858->411), mult. (7332->512), div. (0->0), fcn. (4451->12), ass. (0->206)
t365 = Ifges(5,4) + Ifges(6,4);
t366 = Ifges(5,1) + Ifges(6,1);
t364 = Ifges(6,5) + Ifges(5,5);
t363 = Ifges(5,2) + Ifges(6,2);
t362 = Ifges(6,6) + Ifges(5,6);
t251 = sin(qJ(3));
t255 = cos(qJ(3));
t202 = -mrSges(4,1) * t255 + t251 * mrSges(4,2);
t369 = -mrSges(3,1) + t202;
t250 = sin(qJ(4));
t254 = cos(qJ(4));
t177 = -t250 * t251 + t254 * t255;
t245 = qJD(1) + qJD(2);
t150 = t177 * t245;
t368 = t365 * t150;
t178 = t250 * t255 + t251 * t254;
t151 = t178 * t245;
t367 = t365 * t151;
t334 = t251 / 0.2e1;
t327 = mrSges(5,1) + mrSges(6,1);
t326 = mrSges(5,2) + mrSges(6,2);
t244 = qJD(3) + qJD(4);
t361 = t363 * t150 + t362 * t244 + t367;
t360 = t366 * t151 + t364 * t244 + t368;
t258 = -pkin(8) - pkin(7);
t203 = t258 * t251;
t239 = t255 * pkin(8);
t204 = pkin(7) * t255 + t239;
t130 = t250 * t203 + t254 * t204;
t285 = qJD(3) * t258;
t183 = t251 * t285;
t184 = t255 * t285;
t256 = cos(qJ(2));
t315 = qJD(1) * pkin(1);
t288 = t256 * t315;
t353 = -qJD(4) * t130 + t178 * t288 - t183 * t250 + t254 * t184;
t290 = qJD(4) * t254;
t291 = qJD(4) * t250;
t351 = -t177 * t288 + t254 * t183 + t250 * t184 + t203 * t290 - t204 * t291;
t301 = t245 * t251;
t186 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t301;
t300 = t245 * t255;
t187 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t300;
t359 = (t245 * mrSges(3,2) + t251 * t186 - t255 * t187) * t256;
t310 = qJ(5) * t150;
t252 = sin(qJ(2));
t289 = t252 * t315;
t188 = pkin(7) * t245 + t289;
t281 = pkin(8) * t245 + t188;
t136 = t281 * t255;
t125 = t254 * t136;
t135 = t281 * t251;
t314 = qJD(3) * pkin(3);
t126 = -t135 + t314;
t76 = t126 * t250 + t125;
t42 = t76 + t310;
t276 = t76 * mrSges(5,3) + t42 * mrSges(6,3);
t248 = qJ(3) + qJ(4);
t234 = sin(t248);
t284 = t326 * t234;
t358 = t284 + t369;
t236 = cos(t248);
t357 = t327 * t236;
t356 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t189 = -pkin(2) * t245 - t288;
t275 = mrSges(4,1) * t251 + mrSges(4,2) * t255;
t355 = t189 * t275 + qJD(3) * (Ifges(4,5) * t255 - Ifges(4,6) * t251) / 0.2e1;
t265 = t177 * qJD(4);
t109 = qJD(3) * t177 + t265;
t271 = -qJ(5) * t109 - qJD(5) * t178;
t354 = t271 + t353;
t266 = t178 * qJD(4);
t110 = -qJD(3) * t178 - t266;
t297 = t110 * qJ(5) + t177 * qJD(5);
t352 = t297 + t351;
t142 = t151 * qJ(5);
t123 = t250 * t136;
t75 = t254 * t126 - t123;
t41 = -t142 + t75;
t224 = pkin(1) * t252 + pkin(7);
t325 = -pkin(8) - t224;
t171 = t325 * t251;
t305 = t224 * t255;
t172 = t239 + t305;
t104 = t250 * t171 + t254 * t172;
t350 = t369 * t245;
t348 = t189 * t252 + (t251 ^ 2 + t255 ^ 2) * t188 * t256;
t342 = t150 / 0.2e1;
t340 = t151 / 0.2e1;
t335 = t244 / 0.2e1;
t332 = pkin(1) * t256;
t331 = g(1) * t255;
t249 = qJ(1) + qJ(2);
t235 = sin(t249);
t330 = g(2) * t235;
t323 = mrSges(5,3) * t150;
t322 = mrSges(6,3) * t150;
t321 = Ifges(4,4) * t251;
t320 = Ifges(4,4) * t255;
t317 = Ifges(4,2) * t255;
t316 = pkin(1) * qJD(2);
t286 = qJD(2) * t315;
t308 = qJDD(1) * pkin(1);
t175 = t252 * t308 + t256 * t286;
t242 = qJDD(1) + qJDD(2);
t154 = pkin(7) * t242 + t175;
t292 = qJD(3) * t255;
t99 = -t154 * t251 - t188 * t292;
t312 = t251 * t99;
t293 = qJD(3) * t251;
t98 = t255 * t154 - t188 * t293;
t311 = t255 * t98;
t309 = qJ(5) * t178;
t237 = cos(t249);
t304 = t234 * t237;
t302 = t236 * t237;
t78 = -t254 * t135 - t123;
t223 = pkin(4) * t236;
t226 = pkin(3) * t255 + pkin(2);
t182 = t223 + t226;
t243 = -qJ(5) + t258;
t296 = t235 * t182 + t237 * t243;
t295 = t235 * t226 + t237 * t258;
t294 = t237 * pkin(2) + t235 * pkin(7);
t287 = t256 * t316;
t229 = pkin(3) * t293;
t162 = t242 * t255 - t245 * t293;
t163 = t242 * t251 + t245 * t292;
t68 = t162 * t250 + t163 * t254 + t245 * t265;
t69 = t162 * t254 - t163 * t250 - t245 * t266;
t16 = -t69 * mrSges(6,1) + t68 * mrSges(6,2);
t283 = t326 * t236;
t100 = -pkin(4) * t110 + t229;
t280 = qJD(3) * t325;
t77 = t135 * t250 - t125;
t103 = t254 * t171 - t172 * t250;
t279 = t237 * t182 - t235 * t243;
t129 = t254 * t203 - t204 * t250;
t278 = t237 * t226 - t235 * t258;
t277 = -t302 * t326 - t304 * t327;
t174 = -t252 * t286 + t256 * t308;
t274 = t317 + t321;
t272 = t186 * t255 + t187 * t251;
t147 = -pkin(4) * t177 - t226;
t270 = t284 - t357;
t267 = t251 * (Ifges(4,1) * t255 - t321);
t71 = qJDD(3) * pkin(3) - pkin(8) * t163 + t99;
t74 = pkin(8) * t162 + t98;
t7 = t126 * t290 - t136 * t291 + t250 * t71 + t254 * t74;
t131 = t251 * t280 + t255 * t287;
t132 = -t251 * t287 + t255 * t280;
t24 = t254 * t131 + t250 * t132 + t171 * t290 - t172 * t291;
t153 = -pkin(2) * t242 - t174;
t264 = m(6) * (-pkin(3) * t251 - pkin(4) * t234) - t275;
t152 = -t226 * t245 - t288;
t101 = -pkin(3) * t162 + t153;
t8 = -qJD(4) * t76 - t250 * t74 + t254 * t71;
t25 = -qJD(4) * t104 - t131 * t250 + t254 * t132;
t263 = -qJD(3) * t272 + t255 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t162) - t251 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t163);
t262 = t235 * t356 + t237 * t358 - t302 * t327;
t261 = (m(4) * pkin(7) - t356) * t237 + (-t357 + t358) * t235;
t241 = qJDD(3) + qJDD(4);
t3 = pkin(4) * t241 - qJ(5) * t68 - qJD(5) * t151 + t8;
t36 = pkin(4) * t244 + t41;
t4 = qJ(5) * t69 + qJD(5) * t150 + t7;
t95 = -pkin(4) * t150 + qJD(5) + t152;
t260 = t8 * mrSges(5,1) + t3 * mrSges(6,1) - t7 * mrSges(5,2) - t4 * mrSges(6,2) - t152 * (mrSges(5,1) * t151 + mrSges(5,2) * t150) + t36 * t322 - t95 * (mrSges(6,1) * t151 + mrSges(6,2) * t150) + t75 * t323 + t362 * t69 + t364 * t68 - (t366 * t150 - t367) * t151 / 0.2e1 + t361 * t340 - (t150 * t364 - t151 * t362) * t244 / 0.2e1 + (Ifges(6,3) + Ifges(5,3)) * t241 - (-t363 * t151 + t360 + t368) * t150 / 0.2e1;
t148 = Ifges(4,6) * qJD(3) + t245 * t274;
t208 = Ifges(4,4) * t300;
t149 = Ifges(4,1) * t301 + Ifges(4,5) * qJD(3) + t208;
t23 = -pkin(4) * t69 + qJDD(5) + t101;
t259 = Ifges(3,3) * t242 + t153 * t202 + (t149 + t245 * (-Ifges(4,2) * t251 + t320)) * t292 / 0.2e1 + (t267 * t245 / 0.2e1 + t355) * qJD(3) + t174 * mrSges(3,1) - t175 * mrSges(3,2) + mrSges(4,3) * t311 + t162 * t274 / 0.2e1 + t163 * (Ifges(4,1) * t251 + t320) / 0.2e1 - t148 * t293 / 0.2e1 + (0.2e1 * Ifges(4,5) * t334 + Ifges(4,6) * t255) * qJDD(3) + t255 * (Ifges(4,4) * t163 + Ifges(4,2) * t162) / 0.2e1 + (Ifges(4,1) * t163 + Ifges(4,4) * t162) * t334 + (t361 / 0.2e1 - mrSges(5,1) * t152 - mrSges(6,1) * t95 + t362 * t335 + t365 * t340 + t363 * t342 + t276) * t110 + (-t36 * mrSges(6,3) - t75 * mrSges(5,3) + t360 / 0.2e1 + mrSges(5,2) * t152 + mrSges(6,2) * t95 + t364 * t335 + t366 * t340 + t365 * t342) * t109 + (mrSges(5,2) * t101 + mrSges(6,2) * t23 - mrSges(5,3) * t8 - mrSges(6,3) * t3 + t241 * t364 + t365 * t69 + t366 * t68) * t178 + (-mrSges(5,1) * t101 - mrSges(6,1) * t23 + mrSges(5,3) * t7 + mrSges(6,3) * t4 + t241 * t362 + t363 * t69 + t365 * t68) * t177;
t257 = cos(qJ(1));
t253 = sin(qJ(1));
t240 = t257 * pkin(1);
t238 = t253 * pkin(1);
t230 = t252 * t316;
t227 = -pkin(2) - t332;
t221 = t235 * pkin(2);
t201 = -t226 - t332;
t185 = t230 + t229;
t170 = t177 * qJ(5);
t134 = t147 - t332;
t122 = mrSges(5,1) * t244 - mrSges(5,3) * t151;
t121 = mrSges(6,1) * t244 - mrSges(6,3) * t151;
t120 = -mrSges(5,2) * t244 + t323;
t119 = -mrSges(6,2) * t244 + t322;
t111 = pkin(3) * t301 + pkin(4) * t151;
t102 = -mrSges(4,1) * t162 + mrSges(4,2) * t163;
t97 = t170 + t130;
t96 = t129 - t309;
t93 = -mrSges(5,1) * t150 + mrSges(5,2) * t151;
t92 = -mrSges(6,1) * t150 + mrSges(6,2) * t151;
t89 = t100 + t230;
t80 = t170 + t104;
t79 = t103 - t309;
t55 = -mrSges(5,2) * t241 + mrSges(5,3) * t69;
t54 = -mrSges(6,2) * t241 + mrSges(6,3) * t69;
t53 = mrSges(5,1) * t241 - mrSges(5,3) * t68;
t52 = mrSges(6,1) * t241 - mrSges(6,3) * t68;
t45 = -t142 + t78;
t44 = t77 - t310;
t17 = -mrSges(5,1) * t69 + mrSges(5,2) * t68;
t14 = t25 + t271;
t13 = t24 + t297;
t1 = [t227 * t102 + t201 * t17 + m(5) * (t101 * t201 + t103 * t8 + t104 * t7 + t152 * t185 + t24 * t76 + t25 * t75) + m(6) * (t13 * t42 + t134 * t23 + t14 * t36 + t3 * t79 + t4 * t80 + t89 * t95) + t185 * t93 + t13 * t119 + t24 * t120 + t14 * t121 + t25 * t122 + t134 * t16 + t103 * t53 + t104 * t55 + t89 * t92 + t79 * t52 + t80 * t54 + (-t257 * mrSges(2,1) + t253 * mrSges(2,2) - m(5) * (t240 + t278) - m(4) * (t240 + t294) - m(6) * (t240 + t279) + t262) * g(2) + (-t253 * mrSges(2,1) - t257 * mrSges(2,2) - m(4) * (t221 + t238) - m(5) * (t238 + t295) - m(6) * (t238 + t296) + t261) * g(3) + m(4) * (t153 * t227 - t224 * t312 + t98 * t305) + t263 * t224 + t259 + ((mrSges(3,1) * t256 - mrSges(3,2) * t252) * t242 + (-g(2) * t257 - g(3) * t253 + t174 * t256 + t175 * t252) * m(3) + (m(4) * t348 + t252 * t350 - t359) * qJD(2)) * pkin(1) - mrSges(4,3) * t312 + Ifges(2,3) * qJDD(1); -t226 * t17 + t262 * g(2) + t263 * pkin(7) + t261 * g(3) + (t359 + (-t350 - t92 - t93) * t252) * t315 + t147 * t16 + t129 * t53 + t130 * t55 + t100 * t92 - pkin(2) * t102 + t96 * t52 + t97 * t54 + (-t99 * mrSges(4,3) + t314 * t93) * t251 + t259 + t352 * t119 + t351 * t120 + t354 * t121 + t353 * t122 + (-t279 * g(2) - t296 * g(3) + t147 * t23 + t3 * t96 + t4 * t97 + (t100 - t289) * t95 + t352 * t42 + t354 * t36) * m(6) + (-g(2) * t278 - g(3) * t295 - t101 * t226 + t129 * t8 + t130 * t7 + t351 * t76 + t353 * t75 + (t229 - t289) * t152) * m(5) + (-t294 * g(2) - t221 * g(3) - t348 * t315 - pkin(2) * t153 + (t311 - t312) * pkin(7)) * m(4); (-m(6) * t223 + t202 + t270) * g(1) + (t237 * t264 + t277) * g(3) + Ifges(4,6) * t162 + Ifges(4,5) * t163 - t45 * t119 - t78 * t120 - t44 * t121 - t77 * t122 - t111 * t92 - t98 * mrSges(4,2) + t99 * mrSges(4,1) - m(5) * (t75 * t77 + t76 * t78) + (t234 * t327 - t264 + t283) * t330 + (t148 * t334 + (t317 * t334 - t267 / 0.2e1) * t245 - (t149 + t208) * t255 / 0.2e1 - t355) * t245 - m(6) * (t111 * t95 + t36 * t44 + t42 * t45) + t260 + t272 * t188 + t276 * t151 + (-m(6) * t331 - t93 * t301 + t254 * t53 + (m(6) * t4 + t54 + t55) * t250 + ((m(6) * t42 + t119 + t120) * t254 + (-m(6) * t36 - t121 - t122) * t250) * qJD(4) + (-t331 - t75 * t291 + t76 * t290 + t250 * t7 + t254 * t8 + (-g(3) * t237 - t152 * t245 + t330) * t251) * m(5)) * pkin(3) + Ifges(4,3) * qJDD(3) + (m(6) * t3 + t52) * (pkin(3) * t254 + pkin(4)); -t41 * t119 - t75 * t120 + t42 * t121 + t76 * t122 + t270 * g(1) + pkin(4) * t52 + t260 + (-(-t36 + t41) * t42 + (-g(1) * t236 - g(3) * t304 - t151 * t95 + t3) * pkin(4)) * m(6) + (t283 + (m(6) * pkin(4) + t327) * t234) * t330 + t277 * g(3) + (-pkin(4) * t92 + t276) * t151; -t150 * t119 + t151 * t121 + (g(2) * t237 + g(3) * t235 - t150 * t42 + t151 * t36 + t23) * m(6) + t16;];
tau = t1;
