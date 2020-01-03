% Calculate vector of inverse dynamics joint torques for
% S5RRRPP2
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:30
% EndTime: 2019-12-31 20:51:39
% DurationCPUTime: 4.84s
% Computational Cost: add. (2080->370), mult. (3092->435), div. (0->0), fcn. (1370->8), ass. (0->177)
t324 = Ifges(5,1) + Ifges(6,1);
t322 = Ifges(6,4) + Ifges(5,5);
t323 = Ifges(5,4) + Ifges(4,5);
t314 = -Ifges(6,5) + t323;
t330 = Ifges(6,2) + Ifges(5,3);
t329 = Ifges(4,6) - Ifges(5,6);
t320 = Ifges(5,6) - Ifges(6,6);
t167 = sin(qJ(3));
t170 = cos(qJ(3));
t273 = Ifges(5,5) * t170;
t275 = Ifges(6,4) * t170;
t328 = t324 * t167 - t273 - t275;
t327 = -mrSges(5,1) - mrSges(6,1);
t160 = qJD(1) + qJD(2);
t255 = t160 * t167;
t236 = mrSges(5,2) * t255;
t308 = -mrSges(4,3) * t255 - t236 + (mrSges(4,1) + mrSges(5,1)) * qJD(3);
t254 = t160 * t170;
t243 = mrSges(5,2) * t254;
t110 = qJD(3) * mrSges(5,3) + t243;
t306 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t254 + t110;
t326 = t322 * t255;
t325 = mrSges(6,3) - mrSges(5,2) - mrSges(4,3) + mrSges(3,2);
t319 = t320 * qJD(3) - t330 * t254 + t326;
t318 = -t329 * t167 + t323 * t170;
t168 = sin(qJ(2));
t272 = pkin(1) * qJD(1);
t240 = t168 * t272;
t111 = pkin(7) * t160 + t240;
t248 = qJD(3) * t167;
t171 = cos(qJ(2));
t271 = pkin(1) * qJD(2);
t235 = qJD(1) * t271;
t262 = pkin(1) * qJDD(1);
t101 = t168 * t262 + t171 * t235;
t159 = qJDD(1) + qJDD(2);
t75 = pkin(7) * t159 + t101;
t20 = -t111 * t248 + t170 * t75;
t247 = qJD(3) * t170;
t21 = -t111 * t247 - t167 * t75;
t317 = -t167 * t21 + t170 * t20;
t10 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t20;
t234 = qJDD(4) - t21;
t17 = -qJDD(3) * pkin(3) + t234;
t316 = t10 * t170 + t167 * t17;
t315 = t322 * t167;
t166 = qJ(1) + qJ(2);
t151 = sin(t166);
t152 = cos(t166);
t304 = g(1) * t152 + g(2) * t151;
t130 = Ifges(4,4) * t254;
t313 = Ifges(4,1) * t255 + t314 * qJD(3) + t328 * t160 + t130;
t239 = t171 * t272;
t112 = -pkin(2) * t160 - t239;
t203 = mrSges(4,1) * t167 + mrSges(4,2) * t170;
t153 = t167 * qJ(4);
t221 = pkin(2) + t153;
t290 = pkin(3) + pkin(4);
t25 = t239 + qJD(5) + (t170 * t290 + t221) * t160;
t265 = t170 * mrSges(5,3);
t266 = t170 * mrSges(6,2);
t157 = t170 * pkin(3);
t190 = -t221 - t157;
t38 = t160 * t190 - t239;
t312 = t112 * t203 + t38 * (t167 * mrSges(5,1) - t265) + t25 * (-t167 * mrSges(6,1) + t266);
t94 = t167 * t111;
t67 = -qJD(3) * pkin(3) + qJD(4) + t94;
t163 = qJD(3) * qJ(4);
t95 = t170 * t111;
t77 = t95 + t163;
t193 = -t167 * t77 + t170 * t67;
t88 = t159 * t167 + t160 * t247;
t79 = t88 * mrSges(5,2);
t59 = -qJDD(3) * mrSges(5,1) + t79;
t87 = -t170 * t159 + t160 * t248;
t60 = -mrSges(5,2) * t87 + qJDD(3) * mrSges(5,3);
t311 = t167 * (-qJDD(3) * mrSges(4,1) + mrSges(4,3) * t88 + t59) + t170 * (-qJDD(3) * mrSges(4,2) - mrSges(4,3) * t87 + t60) + m(4) * t317 + m(5) * (qJD(3) * t193 + t316) - t306 * t248 - t308 * t247;
t310 = m(5) + m(6);
t202 = t170 * mrSges(5,1) + t167 * mrSges(5,3);
t83 = t202 * t160;
t154 = t167 * mrSges(6,2);
t251 = t170 * mrSges(6,1) + t154;
t84 = t251 * t160;
t309 = -t83 - t84;
t241 = mrSges(6,3) * t254;
t108 = qJD(3) * mrSges(6,2) - t241;
t307 = t108 + t110;
t49 = qJ(5) * t255 - t94;
t303 = -t49 + qJD(4);
t100 = -t168 * t235 + t171 * t262;
t208 = pkin(2) * t159 + t100;
t302 = m(4) * t208 - mrSges(4,1) * t87 - mrSges(4,2) * t88;
t210 = -m(6) * t290 - mrSges(6,1);
t258 = t151 * t170;
t269 = t167 * mrSges(4,2);
t299 = mrSges(4,1) * t258 + t325 * t152 + (-m(5) * t190 + m(6) * t221 - t170 * t210 + mrSges(3,1) + t154 + t202 - t269) * t151;
t256 = t152 * t170;
t298 = (-mrSges(4,1) + t327) * t256 + t325 * t151 + (-mrSges(3,1) + (mrSges(4,2) - mrSges(6,2) - mrSges(5,3)) * t167) * t152;
t285 = pkin(1) * t171;
t287 = pkin(1) * t168;
t297 = m(4) * pkin(1) * (t112 * t168 + (t167 ^ 2 + t170 ^ 2) * t171 * t111) + (-mrSges(3,1) * t287 - mrSges(3,2) * t285) * t160;
t296 = (-g(1) * t256 - g(2) * t258) * qJ(4);
t169 = sin(qJ(1));
t286 = pkin(1) * t169;
t172 = cos(qJ(1));
t158 = t172 * pkin(1);
t280 = t88 * mrSges(6,3);
t279 = pkin(7) - qJ(5);
t278 = Ifges(4,4) * t167;
t277 = Ifges(4,4) * t170;
t261 = qJ(4) * t170;
t144 = pkin(7) + t287;
t253 = -qJ(5) + t144;
t252 = t152 * pkin(2) + t151 * pkin(7);
t250 = t157 + t153;
t249 = qJD(3) * t160;
t246 = qJD(5) * t167;
t245 = qJD(5) * t170;
t150 = t167 * qJD(4);
t242 = mrSges(6,3) * t255;
t238 = t168 * t271;
t237 = t171 * t271;
t81 = pkin(3) * t248 - qJ(4) * t247 - t150;
t156 = t170 * pkin(4);
t232 = t156 + t250;
t145 = -pkin(2) - t285;
t23 = -t87 * mrSges(6,1) + t88 * mrSges(6,2);
t116 = t279 * t170;
t225 = -t249 / 0.2e1;
t224 = t249 / 0.2e1;
t140 = t152 * pkin(7);
t220 = -pkin(2) * t151 + t140;
t99 = t253 * t170;
t216 = -qJ(5) * t152 + t140;
t211 = pkin(3) * t256 + t152 * t153 + t252;
t50 = -qJ(5) * t254 + t95;
t206 = -qJD(5) + t237;
t115 = -mrSges(4,1) * t170 + t269;
t199 = t170 * Ifges(4,2) + t278;
t196 = Ifges(6,5) * t170 + Ifges(6,6) * t167;
t27 = -qJD(3) * t290 + t303;
t39 = t163 + t50;
t194 = t167 * t27 + t170 * t39;
t96 = t145 - t250;
t186 = t167 * (Ifges(4,1) * t170 - t278);
t185 = t170 * (Ifges(6,2) * t167 + t275);
t184 = t170 * (Ifges(5,3) * t167 + t273);
t183 = (t167 * t67 + t170 * t77) * t171;
t51 = -pkin(4) * t248 - t81;
t181 = pkin(4) * t256 - qJ(5) * t151 + t211;
t178 = qJ(4) * t88 + t150 * t160 + t208;
t2 = -t290 * t87 + qJDD(5) + t178;
t3 = -qJ(5) * t88 - qJDD(3) * t290 - t160 * t246 + t234;
t5 = qJ(5) * t87 - t160 * t245 + t10;
t6 = pkin(3) * t87 - t178;
t70 = Ifges(4,6) * qJD(3) + t160 * t199;
t174 = -t330 * t170 * t87 + (t319 + (t170 * t324 + t315) * t160) * t248 / 0.2e1 + (t160 * (-Ifges(4,2) * t167 + t277) + t313) * t247 / 0.2e1 + (t247 * t67 - t248 * t77 + t316) * mrSges(5,2) + t317 * mrSges(4,3) - t87 * t199 / 0.2e1 - t6 * t202 + (t323 * t167 + t329 * t170) * qJDD(3) / 0.2e1 + ((Ifges(4,1) + t324) * t88 + (-Ifges(4,4) + t322) * t87 + t314 * qJDD(3)) * t167 / 0.2e1 + t170 * (Ifges(4,4) * t88 - Ifges(4,2) * t87 + Ifges(4,6) * qJDD(3)) / 0.2e1 - qJDD(3) * (Ifges(6,5) * t167 - Ifges(6,6) * t170) / 0.2e1 + Ifges(3,3) * t159 + t315 * t87 / 0.2e1 + t100 * mrSges(3,1) - t101 * mrSges(3,2) - (t320 * qJDD(3) + t322 * t88) * t170 / 0.2e1 - t70 * t248 / 0.2e1 + t2 * t251 + (-t167 * t3 - t170 * t5 - t247 * t27 + t39 * t248) * mrSges(6,3) + (t185 + t184) * t225 + t186 * t224 + (t167 * Ifges(4,1) + t277 + t328) * t88 / 0.2e1 - t208 * t115 + (t312 + (t318 / 0.2e1 - t196 / 0.2e1) * qJD(3)) * qJD(3);
t142 = qJ(5) * t248;
t114 = t279 * t167;
t113 = -pkin(2) - t250;
t105 = -qJD(3) * mrSges(6,1) - t242;
t98 = t253 * t167;
t97 = pkin(2) + t232;
t86 = (pkin(3) * t167 - t261) * t160;
t85 = t115 * t160;
t82 = qJD(3) * t116 - t246;
t80 = -pkin(7) * t248 + t142 - t245;
t76 = t156 - t96;
t57 = -qJDD(3) * mrSges(6,1) - t280;
t55 = qJDD(3) * mrSges(6,2) + mrSges(6,3) * t87;
t52 = t81 + t238;
t48 = (-t167 * t290 + t261) * t160;
t31 = t51 - t238;
t29 = qJD(3) * t99 + t167 * t206;
t28 = -t144 * t248 + t170 * t206 + t142;
t22 = mrSges(5,1) * t87 - mrSges(5,3) * t88;
t1 = [m(3) * (t100 * t171 + t101 * t168) * pkin(1) + t85 * t238 + m(6) * (t2 * t76 + t25 * t31 + t27 * t29 + t28 * t39 + t3 * t98 + t5 * t99) + t174 + t96 * t22 + t98 * t57 + t99 * t55 + t29 * t105 + t28 * t108 + t76 * t23 - t52 * t83 + t31 * t84 + Ifges(2,3) * qJDD(1) + m(5) * (t183 * t271 + t38 * t52 + t6 * t96) + (mrSges(3,1) * t285 - mrSges(3,2) * t287) * t159 - t302 * t145 + t297 * qJD(2) + (-m(5) * (t158 + t211) - m(6) * (t158 + t181) - mrSges(2,1) * t172 + mrSges(2,2) * t169 - m(3) * t158 - m(4) * (t158 + t252) + t298) * g(2) + (mrSges(2,1) * t169 + mrSges(2,2) * t172 + m(3) * t286 - m(4) * (t220 - t286) - m(6) * (t216 - t286) - m(5) * (t140 - t286) + t299) * g(1) + (-t167 * t308 + t170 * t306) * t237 + t311 * t144; t82 * t105 + t80 * t108 + t113 * t22 + t114 * t57 + t116 * t55 + t97 * t23 + t51 * t84 - t81 * t83 + t174 + t302 * pkin(2) + (t114 * t3 + t116 * t5 + t2 * t97 + t25 * t51 + t27 * t82 + t39 * t80 - (-t168 * t25 + t171 * t194) * t272) * m(6) + (t113 * t6 + t38 * t81 - (t168 * t38 + t183) * t272) * m(5) + (-t85 - t309) * t240 - t297 * qJD(1) + (-m(4) * t252 - m(5) * t211 - m(6) * t181 + t298) * g(2) + (-m(4) * t220 - m(5) * t140 - m(6) * t216 + t299) * g(1) + ((-t105 + t308) * t167 + (-t108 - t306) * t170) * t239 + t311 * pkin(7); ((-t186 / 0.2e1 + t185 / 0.2e1 + t184 / 0.2e1) * t160 - t312) * t160 + t27 * t241 + t77 * t236 + (Ifges(6,3) + Ifges(5,2) + Ifges(4,3)) * qJDD(3) - (-Ifges(4,2) * t255 + t130 + t313) * t254 / 0.2e1 + t314 * t88 + t318 * t225 + (-t265 - t266 + t203 + (m(5) * pkin(3) + mrSges(5,1) - t210) * t167) * t304 + t306 * t94 + t307 * qJD(4) + t308 * t95 + (-pkin(3) * t17 + qJ(4) * t10 + qJD(4) * t77 - t111 * t193 - t38 * t86 + t296) * m(5) + (-Ifges(4,6) + t320) * t87 + (-m(5) * t250 - m(6) * t232 + t115 - t202 - t251) * g(3) - (t324 * t254 + t319 + t326) * t255 / 0.2e1 - t50 * t105 - t49 * t108 - t48 * t84 + t86 * t83 - pkin(3) * t59 - t20 * mrSges(4,2) + t21 * mrSges(4,1) + t10 * mrSges(5,3) - t17 * mrSges(5,1) - t3 * mrSges(6,1) + t5 * mrSges(6,2) + t70 * t255 / 0.2e1 - t39 * t242 - t67 * t243 + (t55 + t60) * qJ(4) + t196 * t224 + (t5 * qJ(4) - t25 * t48 - t27 * t50 - t290 * t3 + t303 * t39 + t296) * m(6) - t290 * t57; -t280 + t79 + t327 * qJDD(3) - t307 * qJD(3) + t310 * t170 * g(3) + (t309 * t160 - t304 * t310) * t167 + (-qJD(3) * t39 - t25 * t255 + t3) * m(6) + (-qJD(3) * t77 + t255 * t38 + t17) * m(5); (t167 * t105 + t170 * t108) * t160 + (g(1) * t151 - g(2) * t152 + t160 * t194 + t2) * m(6) + t23;];
tau = t1;
