% Calculate vector of inverse dynamics joint torques for
% S5RPRRP11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP11_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:27
% EndTime: 2019-12-31 18:53:51
% DurationCPUTime: 14.04s
% Computational Cost: add. (4408->488), mult. (10598->631), div. (0->0), fcn. (7525->10), ass. (0->205)
t327 = -mrSges(5,3) - mrSges(6,2);
t307 = Ifges(5,1) + Ifges(6,1);
t329 = -Ifges(5,4) + Ifges(6,5);
t306 = Ifges(6,4) + Ifges(5,5);
t305 = Ifges(5,6) - Ifges(6,6);
t321 = Ifges(5,3) + Ifges(6,2);
t148 = sin(pkin(8));
t149 = cos(pkin(8));
t218 = t148 ^ 2 + t149 ^ 2;
t213 = qJD(1) * qJD(2);
t132 = qJDD(1) * qJ(2) + t213;
t151 = sin(qJ(4));
t154 = cos(qJ(4));
t152 = sin(qJ(3));
t155 = cos(qJ(3));
t119 = t148 * t155 + t149 * t152;
t114 = t119 * qJD(1);
t167 = qJD(3) * t154 - t114 * t151;
t118 = t148 * t152 - t149 * t155;
t115 = t118 * qJD(3);
t84 = -qJD(1) * t115 + qJDD(1) * t119;
t53 = qJD(4) * t167 + qJDD(3) * t151 + t154 * t84;
t279 = t53 / 0.2e1;
t94 = qJD(3) * t151 + t114 * t154;
t54 = qJD(4) * t94 - qJDD(3) * t154 + t151 * t84;
t277 = t54 / 0.2e1;
t116 = t119 * qJD(3);
t85 = -qJD(1) * t116 - qJDD(1) * t118;
t78 = qJDD(4) - t85;
t276 = t78 / 0.2e1;
t310 = -m(5) - m(6);
t328 = -m(4) + t310;
t187 = -mrSges(3,1) * t149 + mrSges(3,2) * t148;
t326 = -m(3) * pkin(1) - mrSges(2,1) + t187;
t182 = t154 * mrSges(6,1) + t151 * mrSges(6,3);
t184 = mrSges(5,1) * t154 - mrSges(5,2) * t151;
t325 = -t182 - t184;
t243 = mrSges(4,3) * t114;
t299 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t167 - mrSges(5,2) * t94 - t243;
t248 = pkin(6) + qJ(2);
t125 = t248 * t148;
t120 = qJD(1) * t125;
t126 = t248 * t149;
t121 = qJD(1) * t126;
t86 = -t120 * t155 - t121 * t152;
t79 = -qJD(3) * pkin(3) - t86;
t324 = -m(5) * t79 + t299;
t113 = t118 * qJD(1);
t318 = qJD(4) + t113;
t322 = t306 * t78 + t307 * t53 + t329 * t54;
t309 = mrSges(5,1) + mrSges(6,1);
t308 = mrSges(5,2) - mrSges(6,3);
t320 = t167 * t305 + t306 * t94 + t321 * t318;
t261 = Ifges(6,5) * t167;
t92 = Ifges(5,4) * t167;
t302 = t306 * t318 + t307 * t94 - t261 + t92;
t292 = -t151 * t305 + t154 * t306;
t239 = Ifges(6,5) * t151;
t241 = Ifges(5,4) * t151;
t291 = t154 * t307 + t239 - t241;
t214 = qJD(4) * t154;
t162 = -t115 * t151 + t119 * t214;
t319 = -t305 * t54 + t306 * t53 + t321 * t78;
t317 = Ifges(6,5) * t279 + Ifges(6,6) * t276 - t53 * Ifges(5,4) / 0.2e1 - t78 * Ifges(5,6) / 0.2e1 + (Ifges(6,3) + Ifges(5,2)) * t277;
t196 = m(3) * qJ(2) + mrSges(3,3);
t316 = -t196 - mrSges(4,3) + mrSges(2,2);
t147 = pkin(8) + qJ(3);
t141 = sin(t147);
t142 = cos(t147);
t186 = t142 * mrSges(4,1) - t141 * mrSges(4,2);
t315 = t141 * t327 - t186;
t313 = m(6) * pkin(4) + t309;
t312 = m(6) * qJ(5) - t308;
t215 = qJD(4) * t151;
t137 = pkin(2) * t149 + pkin(1);
t122 = -qJDD(1) * t137 + qJDD(2);
t32 = -pkin(3) * t85 - pkin(7) * t84 + t122;
t191 = pkin(6) * qJDD(1) + t132;
t100 = t191 * t148;
t101 = t191 * t149;
t216 = qJD(3) * t155;
t217 = qJD(3) * t152;
t35 = -t100 * t152 + t101 * t155 - t120 * t216 - t121 * t217;
t33 = qJDD(3) * pkin(7) + t35;
t123 = -qJD(1) * t137 + qJD(2);
t65 = pkin(3) * t113 - pkin(7) * t114 + t123;
t87 = -t120 * t152 + t121 * t155;
t80 = qJD(3) * pkin(7) + t87;
t3 = t151 * t32 + t154 * t33 + t214 * t65 - t215 * t80;
t1 = qJ(5) * t78 + qJD(5) * t318 + t3;
t26 = t151 * t65 + t154 * t80;
t4 = -qJD(4) * t26 - t151 * t33 + t154 * t32;
t2 = -pkin(4) * t78 + qJDD(5) - t4;
t311 = t4 * mrSges(5,1) - t2 * mrSges(6,1) - t3 * mrSges(5,2) + t1 * mrSges(6,3);
t17 = -mrSges(6,2) * t54 + mrSges(6,3) * t78;
t20 = -mrSges(5,2) * t78 - mrSges(5,3) * t54;
t304 = t17 + t20;
t18 = mrSges(5,1) * t78 - mrSges(5,3) * t53;
t19 = -mrSges(6,1) * t78 + mrSges(6,2) * t53;
t303 = t19 - t18;
t61 = mrSges(6,2) * t167 + mrSges(6,3) * t318;
t62 = -mrSges(5,2) * t318 + mrSges(5,3) * t167;
t301 = t61 + t62;
t63 = mrSges(5,1) * t318 - mrSges(5,3) * t94;
t64 = -mrSges(6,1) * t318 + mrSges(6,2) * t94;
t300 = t63 - t64;
t173 = pkin(4) * t151 - qJ(5) * t154;
t298 = -qJD(5) * t151 + t173 * t318 - t87;
t297 = Ifges(4,5) * qJD(3);
t296 = Ifges(4,6) * qJD(3);
t295 = -t125 * t155 - t126 * t152;
t294 = pkin(3) * t142 + pkin(7) * t141;
t181 = mrSges(6,1) * t151 - mrSges(6,3) * t154;
t183 = mrSges(5,1) * t151 + mrSges(5,2) * t154;
t27 = -pkin(4) * t167 - qJ(5) * t94 + t79;
t293 = t181 * t27 + t183 * t79;
t25 = -t151 * t80 + t154 * t65;
t290 = -t25 + qJD(5);
t22 = qJ(5) * t318 + t26;
t289 = -t22 * mrSges(6,2) - t26 * mrSges(5,3);
t21 = -pkin(4) * t318 + t290;
t288 = t21 * mrSges(6,2) - t25 * mrSges(5,3);
t287 = -t151 * t4 + t154 * t3;
t286 = t1 * t154 + t151 * t2;
t282 = pkin(7) * t310;
t174 = pkin(4) * t154 + qJ(5) * t151;
t124 = -pkin(3) - t174;
t281 = (mrSges(4,2) + t327) * t142 + (m(5) * pkin(3) - m(6) * t124 + mrSges(4,1) - t325) * t141;
t83 = pkin(3) * t118 - pkin(7) * t119 - t137;
t90 = -t125 * t152 + t126 * t155;
t245 = t151 * t83 + t154 * t90;
t66 = -qJD(2) * t118 + qJD(3) * t295;
t82 = pkin(3) * t116 + pkin(7) * t115;
t13 = -qJD(4) * t245 - t151 * t66 + t154 * t82;
t278 = -t54 / 0.2e1;
t275 = t167 / 0.2e1;
t274 = -t167 / 0.2e1;
t273 = -t94 / 0.2e1;
t272 = t94 / 0.2e1;
t271 = -t318 / 0.2e1;
t269 = t113 / 0.2e1;
t267 = t114 / 0.2e1;
t264 = t151 / 0.2e1;
t262 = Ifges(5,4) * t94;
t258 = g(3) * t141;
t81 = pkin(3) * t114 + pkin(7) * t113;
t38 = t151 * t81 + t154 * t86;
t244 = mrSges(4,3) * t113;
t242 = Ifges(4,4) * t114;
t240 = Ifges(5,4) * t154;
t238 = Ifges(6,5) * t154;
t233 = t113 * t151;
t232 = t113 * t154;
t156 = cos(qJ(1));
t227 = t141 * t156;
t226 = t142 * t156;
t153 = sin(qJ(1));
t224 = t151 * t153;
t223 = t153 * t154;
t222 = t154 * t115;
t221 = t154 * t156;
t220 = t156 * t151;
t211 = qJDD(1) * t148;
t210 = qJDD(1) * t149;
t91 = Ifges(6,5) * t94;
t40 = Ifges(6,6) * t318 - Ifges(6,3) * t167 + t91;
t204 = t40 * t264;
t195 = -mrSges(4,1) * t85 + mrSges(4,2) * t84;
t193 = -t215 / 0.2e1;
t192 = t214 / 0.2e1;
t190 = t137 * t156 + t153 * t248;
t188 = -mrSges(3,1) * t210 + mrSges(3,2) * t211;
t178 = -Ifges(5,2) * t151 + t240;
t175 = Ifges(6,3) * t151 + t238;
t37 = -t151 * t86 + t154 * t81;
t51 = -t151 * t90 + t154 * t83;
t36 = -t100 * t155 - t101 * t152 + t120 * t217 - t121 * t216;
t12 = t151 * t82 + t154 * t66 + t214 * t83 - t215 * t90;
t161 = t119 * t215 + t222;
t34 = -qJDD(3) * pkin(3) - t36;
t67 = qJD(2) * t119 + qJD(3) * t90;
t140 = -qJDD(1) * pkin(1) + qJDD(2);
t110 = t142 * t221 + t224;
t109 = t142 * t220 - t223;
t108 = t142 * t223 - t220;
t107 = t142 * t224 + t221;
t102 = Ifges(4,4) * t113;
t95 = -qJD(3) * mrSges(4,2) - t244;
t70 = t114 * Ifges(4,1) - t102 + t297;
t69 = -t113 * Ifges(4,2) + t242 + t296;
t57 = -mrSges(6,1) * t167 - mrSges(6,3) * t94;
t56 = pkin(4) * t94 - qJ(5) * t167;
t55 = t119 * t173 - t295;
t43 = Ifges(5,2) * t167 + Ifges(5,6) * t318 + t262;
t30 = -pkin(4) * t118 - t51;
t29 = qJ(5) * t118 + t245;
t24 = -pkin(4) * t114 - t37;
t23 = qJ(5) * t114 + t38;
t16 = -t173 * t115 + (qJD(4) * t174 - qJD(5) * t154) * t119 + t67;
t15 = mrSges(5,1) * t54 + mrSges(5,2) * t53;
t14 = mrSges(6,1) * t54 - mrSges(6,3) * t53;
t11 = -pkin(4) * t116 - t13;
t6 = qJ(5) * t116 + qJD(5) * t118 + t12;
t5 = pkin(4) * t54 - qJ(5) * t53 - qJD(5) * t94 + t34;
t7 = [(-Ifges(4,1) * t115 - Ifges(4,4) * t116) * t267 - t113 * (-Ifges(4,4) * t115 - Ifges(4,2) * t116) / 0.2e1 + qJD(3) * (-Ifges(4,5) * t115 - Ifges(4,6) * t116) / 0.2e1 + t123 * (mrSges(4,1) * t116 - mrSges(4,2) * t115) - t162 * t43 / 0.2e1 + t245 * t20 + 0.2e1 * t218 * t132 * mrSges(3,3) + (-m(4) * t190 + t327 * t227 + t310 * (pkin(3) * t226 + pkin(7) * t227 + t190) - t313 * t110 - t312 * t109 + (-t186 + t326) * t156 + t316 * t153) * g(2) + (t313 * t108 + t312 * t107 + (m(4) * t137 + t310 * (-t137 - t294) - t315 - t326) * t153 + (t248 * t328 + t316) * t156) * g(1) + (-m(4) * t86 - t324) * t67 + (t319 / 0.2e1 + t122 * mrSges(4,1) - t35 * mrSges(4,3) - Ifges(4,4) * t84 - Ifges(4,2) * t85 - Ifges(4,6) * qJDD(3) + Ifges(5,6) * t278 + Ifges(6,6) * t277 + t321 * t276 + t306 * t279 + t311) * t118 + (Ifges(3,4) * t148 + Ifges(3,2) * t149) * t210 + (Ifges(3,1) * t148 + Ifges(3,4) * t149) * t211 + (-Ifges(6,5) * t161 + Ifges(6,6) * t116 + Ifges(6,3) * t162) * t274 + (-Ifges(5,4) * t161 - Ifges(5,2) * t162 + Ifges(5,6) * t116) * t275 + m(3) * (-pkin(1) * t140 + (t132 + t213) * qJ(2) * t218) - t137 * t195 + (t306 * t116 - t307 * t161 + t162 * t329) * t272 - t302 * t222 / 0.2e1 + t320 * t116 / 0.2e1 + m(6) * (t1 * t29 + t11 * t21 + t16 * t27 + t2 * t30 + t22 * t6 + t5 * t55) + t66 * t95 + t16 * t57 + t6 * t61 + t12 * t62 + t13 * t63 + t11 * t64 + t51 * t18 + t55 * t14 + t29 * t17 + t30 * t19 + (t115 * t86 - t116 * t87 - t295 * t84 + t85 * t90) * mrSges(4,3) + (mrSges(4,1) * t295 - mrSges(4,2) * t90) * qJDD(3) + m(5) * (t12 * t26 + t13 * t25 + t245 * t3 - t295 * t34 + t4 * t51) + m(4) * (-t122 * t137 + t295 * t36 + t35 * t90 + t66 * t87) - t295 * t15 + (t116 * t321 - t161 * t306 - t162 * t305) * t318 / 0.2e1 + (t122 * mrSges(4,2) - t36 * mrSges(4,3) + Ifges(4,1) * t84 + Ifges(4,4) * t85 + Ifges(4,5) * qJDD(3) + t175 * t277 + t178 * t278 + t181 * t5 + t183 * t34 + t192 * t40 + t292 * t276 + t291 * t279 + (-mrSges(6,2) * t1 - mrSges(5,3) * t3 + t317) * t151 + t302 * t193 + (t322 / 0.2e1 + t2 * mrSges(6,2) - t4 * mrSges(5,3)) * t154) * t119 + t25 * (mrSges(5,1) * t116 + mrSges(5,3) * t161) + t21 * (-mrSges(6,1) * t116 - mrSges(6,2) * t161) + t140 * t187 - pkin(1) * t188 + t27 * (mrSges(6,1) * t162 + mrSges(6,3) * t161) + t79 * (mrSges(5,1) * t162 - mrSges(5,2) * t161) + t22 * (-mrSges(6,2) * t162 + mrSges(6,3) * t116) + t26 * (-mrSges(5,2) * t116 - mrSges(5,3) * t162) - t115 * t70 / 0.2e1 - t116 * t69 / 0.2e1 - t115 * t204 + Ifges(2,3) * qJDD(1); t113 * t95 + (-t57 + t299) * t114 + (t301 * t318 - t303) * t154 + (-t300 * t318 + t304) * t151 + m(3) * t140 + t188 + t195 + (-g(1) * t153 + g(2) * t156) * (m(3) - t328) - t196 * t218 * qJD(1) ^ 2 + (t1 * t151 - t114 * t27 - t2 * t154 + t318 * (t151 * t21 + t154 * t22)) * m(6) + (-t114 * t79 + t3 * t151 + t4 * t154 + t318 * (-t151 * t25 + t154 * t26)) * m(5) + (t113 * t87 + t114 * t86 + t122) * m(4); t322 * t264 + (-t238 + t240) * t279 + t239 * t277 + t241 * t278 + (t192 + t232 / 0.2e1) * t302 + (-t244 - t95) * t86 + (t243 + t324) * t87 + (t310 * t294 + (-m(6) * t174 + t325) * t142 + t315) * g(3) + (Ifges(5,2) * t278 - Ifges(6,3) * t277 + t305 * t276 - t317) * t154 + t69 * t267 + (t21 * t232 - t22 * t233 + t286) * mrSges(6,2) + (-t232 * t25 - t233 * t26 + t287) * mrSges(5,3) + t288 * t214 + t289 * t215 + (t156 * t281 + t226 * t282) * g(1) + (-pkin(3) * t34 - t25 * t37 - t26 * t38) * m(5) + (t124 * t5 - t21 * t24 - t22 * t23 + t298 * t27) * m(6) + (((-t151 * t26 - t154 * t25) * qJD(4) + t287) * m(5) + ((-t151 * t22 + t154 * t21) * qJD(4) + t286) * m(6) - t300 * t214 - t301 * t215 + t303 * t151 + t304 * t154) * pkin(7) + (t204 + t293 - (-t178 / 0.2e1 + t175 / 0.2e1) * t167) * qJD(4) + (-t102 + t70) * t269 - (-Ifges(4,1) * t113 - t242 + t320) * t114 / 0.2e1 + (-Ifges(4,2) * t269 + Ifges(5,6) * t274 + Ifges(6,6) * t275 - t22 * mrSges(6,3) - t25 * mrSges(5,1) + t21 * mrSges(6,1) + t26 * mrSges(5,2) + t296 / 0.2e1 - t123 * mrSges(4,1) + t306 * t273 + t321 * t271) * t114 - (t178 * t274 + t175 * t275 - t297 / 0.2e1 - t123 * mrSges(4,2) + t291 * t273 + t292 * t271 - t293) * t113 + t298 * t57 + (t142 * t282 + t281) * g(2) * t153 + Ifges(4,6) * t85 + Ifges(4,5) * t84 - t23 * t61 - t38 * t62 - t37 * t63 - t24 * t64 - t35 * mrSges(4,2) + t36 * mrSges(4,1) - pkin(3) * t15 + (t291 * t94 + t292 * t318) * qJD(4) / 0.2e1 + (t276 * t306 + t279 * t307) * t151 - (t43 / 0.2e1 - t40 / 0.2e1) * t233 - t34 * t184 - t5 * t182 + t124 * t14 + t43 * t193 + Ifges(4,3) * qJDD(3); (t183 + t181) * t258 + t43 * t272 + (-t79 * mrSges(5,1) - t27 * mrSges(6,1) + Ifges(6,3) * t275 - t289) * t94 + (t173 * t258 - t21 * t26 - t27 * t56 - pkin(4) * t2 + qJ(5) * t1 - g(2) * (-pkin(4) * t107 + qJ(5) * t108) - g(1) * (-pkin(4) * t109 + qJ(5) * t110) + t290 * t22) * m(6) + t261 * t275 + t311 - (t79 * mrSges(5,2) - t27 * mrSges(6,3) + t288) * t167 + (t167 * t306 - t305 * t94) * t271 + (t167 * t307 - t262 + t40 + t91) * t273 + (-Ifges(5,2) * t94 + t302 + t92) * t274 + t300 * t26 - t301 * t25 - t56 * t57 + qJD(5) * t61 + qJ(5) * t17 - pkin(4) * t19 + (t107 * t309 + t108 * t308) * g(2) + (t109 * t309 + t110 * t308) * g(1) + t319; -t318 * t61 + t94 * t57 + (-g(1) * t109 - g(2) * t107 - t151 * t258 - t22 * t318 + t27 * t94 + t2) * m(6) + t19;];
tau = t7;
