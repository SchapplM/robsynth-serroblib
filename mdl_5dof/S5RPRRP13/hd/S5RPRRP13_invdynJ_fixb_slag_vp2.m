% Calculate vector of inverse dynamics joint torques for
% S5RPRRP13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP13_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP13_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP13_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP13_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP13_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:58:32
% EndTime: 2019-12-31 18:58:50
% DurationCPUTime: 10.37s
% Computational Cost: add. (2653->462), mult. (5237->605), div. (0->0), fcn. (2833->6), ass. (0->217)
t332 = Ifges(5,1) + Ifges(6,1);
t317 = Ifges(5,5) + Ifges(6,4);
t331 = Ifges(5,6) - Ifges(6,6);
t123 = sin(qJ(4));
t126 = cos(qJ(4));
t216 = t126 * qJD(3);
t127 = cos(qJ(3));
t225 = qJD(1) * t127;
t88 = t123 * t225 - t216;
t242 = qJD(4) * t88;
t124 = sin(qJ(3));
t214 = qJD(1) * qJD(3);
t95 = qJDD(1) * t127 - t124 * t214;
t39 = qJDD(3) * t123 + t126 * t95 - t242;
t285 = t39 / 0.2e1;
t224 = qJD(3) * t123;
t89 = t126 * t225 + t224;
t40 = qJD(4) * t89 - t126 * qJDD(3) + t123 * t95;
t283 = t40 / 0.2e1;
t96 = -qJDD(1) * t124 - t127 * t214;
t86 = qJDD(4) - t96;
t282 = t86 / 0.2e1;
t280 = t88 / 0.2e1;
t334 = -t89 / 0.2e1;
t226 = qJD(1) * t124;
t108 = qJD(4) + t226;
t333 = -t108 / 0.2e1;
t330 = Ifges(5,3) + Ifges(6,2);
t305 = -t331 * t123 + t317 * t126;
t251 = Ifges(6,5) * t123;
t253 = Ifges(5,4) * t123;
t303 = t332 * t126 + t251 - t253;
t284 = -t40 / 0.2e1;
t329 = t317 * t282 + (-Ifges(5,4) + Ifges(6,5)) * t283 + t332 * t285;
t171 = t126 * mrSges(6,1) + t123 * mrSges(6,3);
t173 = mrSges(5,1) * t126 - mrSges(5,2) * t123;
t328 = -t171 - t173;
t219 = qJD(4) * t126;
t220 = qJD(4) * t123;
t215 = qJD(1) * qJD(2);
t107 = qJDD(1) * qJ(2) + t215;
t42 = -pkin(3) * t96 - pkin(7) * t95 + t107;
t129 = -pkin(1) - pkin(6);
t104 = qJDD(1) * t129 + qJDD(2);
t106 = qJD(1) * t129 + qJD(2);
t222 = qJD(3) * t127;
t55 = t124 * t104 + t106 * t222;
t51 = qJDD(3) * pkin(7) + t55;
t120 = t127 * pkin(7);
t98 = pkin(3) * t124 + qJ(2) - t120;
t70 = t98 * qJD(1);
t97 = t124 * t106;
t76 = qJD(3) * pkin(7) + t97;
t4 = t123 * t42 + t126 * t51 + t70 * t219 - t220 * t76;
t31 = t123 * t70 + t126 * t76;
t5 = -qJD(4) * t31 - t123 * t51 + t126 * t42;
t176 = -t5 * t123 + t4 * t126;
t30 = -t123 * t76 + t126 * t70;
t327 = -t30 * t219 - t31 * t220 + t176;
t1 = qJ(5) * t86 + qJD(5) * t108 + t4;
t2 = -pkin(4) * t86 + qJDD(5) - t5;
t177 = t1 * t126 + t123 * t2;
t21 = -pkin(4) * t108 + qJD(5) - t30;
t22 = qJ(5) * t108 + t31;
t326 = t21 * t219 - t22 * t220 + t177;
t14 = mrSges(5,1) * t86 - mrSges(5,3) * t39;
t15 = -t86 * mrSges(6,1) + t39 * mrSges(6,2);
t16 = -mrSges(5,2) * t86 - mrSges(5,3) * t40;
t17 = -mrSges(6,2) * t40 + mrSges(6,3) * t86;
t325 = (t16 + t17) * t126 + (-t14 + t15) * t123;
t255 = Ifges(4,4) * t124;
t169 = t127 * Ifges(4,1) - t255;
t82 = Ifges(6,5) * t89;
t24 = t108 * Ifges(6,6) + t88 * Ifges(6,3) + t82;
t269 = Ifges(6,5) * t88;
t83 = Ifges(5,4) * t88;
t314 = t317 * t108 + t332 * t89 + t269 - t83;
t324 = Ifges(4,5) * qJD(3) + qJD(1) * t169 + t123 * t24 + t126 * t314;
t254 = Ifges(4,4) * t127;
t164 = -t124 * Ifges(4,2) + t254;
t323 = t317 * t334 + t331 * t280 + t330 * t333 + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t164 / 0.2e1;
t322 = t95 / 0.2e1;
t321 = t96 / 0.2e1;
t320 = m(5) + m(6);
t318 = mrSges(5,3) + mrSges(6,2);
t11 = mrSges(5,1) * t40 + mrSges(5,2) * t39;
t316 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t95 - t11;
t271 = mrSges(5,3) * t88;
t56 = -mrSges(5,2) * t108 - t271;
t273 = mrSges(6,2) * t88;
t59 = mrSges(6,3) * t108 - t273;
t257 = t56 + t59;
t270 = mrSges(5,3) * t89;
t57 = mrSges(5,1) * t108 - t270;
t272 = mrSges(6,2) * t89;
t58 = -mrSges(6,1) * t108 + t272;
t256 = -t57 + t58;
t153 = pkin(4) * t123 - qJ(5) * t126;
t313 = -qJD(5) * t123 + t108 * t153 - t97;
t235 = t124 * t129;
t312 = t123 * t98 + t126 * t235;
t311 = qJD(3) * mrSges(4,1) - mrSges(5,1) * t88 - mrSges(5,2) * t89 - mrSges(4,3) * t225;
t125 = sin(qJ(1));
t238 = t124 * t125;
t233 = t125 * t127;
t154 = pkin(4) * t126 + qJ(5) * t123;
t99 = -pkin(3) - t154;
t309 = m(5) * pkin(3) - m(6) * t99 - t328;
t308 = -t305 * t124 + t330 * t127;
t307 = -t303 * t124 + t317 * t127;
t306 = t317 * t123 + t331 * t126;
t250 = Ifges(6,5) * t126;
t252 = Ifges(5,4) * t126;
t304 = t332 * t123 - t250 + t252;
t300 = t317 * t39 + t330 * t86 - t331 * t40;
t223 = qJD(3) * t124;
t54 = t104 * t127 - t106 * t223;
t299 = t124 * t55 + t127 * t54;
t175 = mrSges(4,1) * t127 - mrSges(4,2) * t124;
t298 = t127 * (-Ifges(4,1) * t124 - t254) / 0.2e1 + qJ(2) * t175;
t297 = -t39 * Ifges(6,5) / 0.2e1 - t86 * Ifges(6,6) / 0.2e1 + Ifges(5,4) * t285 + Ifges(5,6) * t282 + (Ifges(6,3) + Ifges(5,2)) * t284;
t296 = -mrSges(4,3) - mrSges(2,1) + mrSges(3,2);
t295 = g(1) * t238;
t174 = mrSges(4,1) * t124 + mrSges(4,2) * t127;
t294 = -t174 + mrSges(2,2) - mrSges(3,3);
t179 = pkin(3) * t127 + pkin(7) * t124;
t87 = qJD(3) * t179 + qJD(2);
t290 = -qJD(4) * t312 + t126 * t87;
t289 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t288 = m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3);
t286 = t5 * mrSges(5,1) - t2 * mrSges(6,1) - t4 * mrSges(5,2) + t1 * mrSges(6,3);
t130 = qJD(1) ^ 2;
t281 = -t88 / 0.2e1;
t278 = t89 / 0.2e1;
t266 = g(1) * t125;
t128 = cos(qJ(1));
t265 = g(2) * t128;
t260 = t89 * Ifges(5,4);
t259 = -qJD(1) / 0.2e1;
t258 = qJD(4) / 0.2e1;
t232 = t126 * t127;
t93 = t179 * qJD(1);
t45 = t106 * t232 + t123 * t93;
t246 = t126 * t93;
t245 = t126 * t98;
t241 = t106 * t127;
t240 = t123 * t124;
t239 = t123 * t127;
t237 = t124 * t126;
t236 = t124 * t128;
t234 = t125 * t126;
t231 = t127 * t128;
t229 = t128 * t126;
t228 = pkin(3) * t233 + pkin(7) * t238;
t227 = t128 * pkin(1) + t125 * qJ(2);
t221 = qJD(3) * t129;
t218 = qJD(4) * t127;
t217 = qJDD(1) * mrSges(3,2);
t199 = t127 * t221;
t210 = t123 * t87 + t126 * t199 + t98 * t219;
t202 = t128 * pkin(6) + t227;
t201 = t123 * t223;
t198 = t129 * t220;
t188 = t219 / 0.2e1;
t187 = -t218 / 0.2e1;
t185 = t123 * t129 - pkin(4);
t184 = -t214 / 0.2e1;
t183 = (t107 + t215) * qJ(2);
t178 = t265 - t266;
t172 = mrSges(5,1) * t123 + mrSges(5,2) * t126;
t170 = t123 * mrSges(6,1) - t126 * mrSges(6,3);
t163 = -Ifges(5,2) * t123 + t252;
t162 = Ifges(5,2) * t126 + t253;
t159 = -Ifges(4,5) * t124 - Ifges(4,6) * t127;
t156 = Ifges(6,3) * t123 + t250;
t155 = -Ifges(6,3) * t126 + t251;
t152 = t22 * t123 - t21 * t126;
t151 = t31 * t123 + t30 * t126;
t77 = -qJD(3) * pkin(3) - t241;
t149 = -t129 + t153;
t50 = -qJDD(3) * pkin(3) - t54;
t146 = t124 * (-Ifges(4,2) * t127 - t255);
t143 = -t130 * qJ(2) + t178;
t141 = t123 * t218 + t124 * t216;
t140 = -t126 * t218 + t201;
t139 = -t123 * t257 + t126 * t256;
t134 = Ifges(5,6) * t127 - t124 * t163;
t133 = Ifges(6,6) * t127 - t124 * t156;
t119 = t128 * qJ(2);
t114 = -qJDD(1) * pkin(1) + qJDD(2);
t101 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t226;
t92 = t174 * qJD(1);
t81 = t172 * t127;
t74 = -t123 * t125 + t124 * t229;
t73 = t123 * t236 + t234;
t72 = t123 * t128 + t124 * t234;
t71 = t123 * t238 - t229;
t65 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t96;
t62 = t149 * t127;
t60 = -t123 * t235 + t245;
t53 = t124 * t185 - t245;
t52 = qJ(5) * t124 + t312;
t47 = mrSges(6,1) * t88 - mrSges(6,3) * t89;
t46 = pkin(4) * t89 + qJ(5) * t88;
t44 = -t106 * t239 + t246;
t41 = -t246 + (-pkin(4) * qJD(1) + t106 * t123) * t127;
t38 = qJ(5) * t225 + t45;
t27 = -t88 * Ifges(5,2) + t108 * Ifges(5,6) + t260;
t23 = pkin(4) * t88 - qJ(5) * t89 + t77;
t20 = (qJD(4) * t154 - qJD(5) * t126) * t127 - t149 * t223;
t19 = -t123 * t199 + t290;
t18 = -t124 * t198 + t210;
t13 = t185 * t222 - t290;
t12 = qJ(5) * t222 + (qJD(5) - t198) * t124 + t210;
t10 = mrSges(6,1) * t40 - mrSges(6,3) * t39;
t3 = pkin(4) * t40 - qJ(5) * t39 - qJD(5) * t89 + t50;
t6 = [(-Ifges(4,4) * t95 / 0.2e1 - Ifges(4,2) * t96 / 0.2e1 + Ifges(6,6) * t283 + Ifges(5,6) * t284 + t317 * t285 + t330 * t282 + t286 + t300 / 0.2e1 - Ifges(4,6) * qJDD(3) - t311 * t221) * t124 + (qJD(3) * t307 - t218 * t304) * t278 + (qJD(3) * t308 - t218 * t306) * t108 / 0.2e1 + t232 * t329 + t298 * t214 + m(4) * (t129 * t299 + t183) - t299 * mrSges(4,3) + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) + t114 * mrSges(3,2) + qJD(2) * t92 + qJ(2) * (-mrSges(4,1) * t96 + mrSges(4,2) * t95) + t50 * t81 + t23 * (-mrSges(6,1) * t140 + mrSges(6,3) * t141) + t77 * (-mrSges(5,1) * t140 - mrSges(5,2) * t141) - t297 * t239 + t62 * t10 + t52 * t17 + t53 * t15 + t18 * t56 + t19 * t57 + t13 * t58 + t12 * t59 + t60 * t14 + t20 * t47 + m(5) * (t129 * t223 * t77 + t18 * t31 + t19 * t30 + t312 * t4 + t5 * t60) + t312 * t16 + (t174 + 0.2e1 * mrSges(3,3)) * t107 + (-t1 * t239 + t140 * t22 - t141 * t21 + t2 * t232) * mrSges(6,2) + (Ifges(4,1) * t322 + Ifges(4,4) * t321 + Ifges(4,5) * qJDD(3) + t156 * t283 + t163 * t284 + t3 * t170 + t24 * t188 + t305 * t282 + t303 * t285 + (-m(5) * t50 + t316) * t129) * t127 + m(6) * (t1 * t52 + t12 * t22 + t13 * t21 + t2 * t53 + t20 * t23 + t3 * t62) + (t201 / 0.2e1 + t126 * t187) * t27 + (-m(3) * t227 - m(4) * t202 - t320 * (pkin(3) * t238 - pkin(7) * t233 + t202) - t289 * t72 - t288 * t71 + t318 * t233 + t296 * t128 + t294 * t125) * g(2) + (-t320 * (pkin(3) * t236 - pkin(7) * t231 + t125 * t129 + t119) - t289 * t74 - t288 * t73 + t318 * t231 + (-m(3) - m(4)) * t119 + t294 * t128 + (m(3) * pkin(1) - m(4) * t129 - t296) * t125) * g(1) + m(3) * (-pkin(1) * t114 + t183) - t324 * t223 / 0.2e1 + (t30 * mrSges(5,1) - t21 * mrSges(6,1) - t31 * mrSges(5,2) + t22 * mrSges(6,3) - t323) * t222 + t314 * t123 * t187 + t101 * t199 - pkin(1) * t217 + t146 * t184 + (t31 * t140 + t30 * t141 - t5 * t232 - t4 * t239) * mrSges(5,3) + t65 * t235 + t164 * t321 + t169 * t322 + (qJD(3) * t133 - t155 * t218) * t280 + (qJD(3) * t134 - t162 * t218) * t281 + qJD(3) ^ 2 * t159 / 0.2e1; t217 - t130 * mrSges(3,3) + t143 * m(4) + (t114 + t143) * m(3) + (-m(5) * t151 - m(6) * t152 + t139 - t92) * qJD(1) + (-t10 + (t123 * t256 + t126 * t257 + t101) * qJD(3) + m(5) * (t216 * t31 - t224 * t30 - t50) + m(6) * (t21 * t224 + t216 * t22 - t3) + m(4) * t54 + t316) * t127 + (t65 + (t47 - t311) * qJD(3) + t139 * qJD(4) + m(5) * (qJD(3) * t77 + t327) + m(6) * (qJD(3) * t23 + t326) + m(4) * t55 + t325) * t124 + t320 * t178; t314 * t188 + t304 * t285 + t306 * t282 + t123 * t329 + t159 * t184 + t297 * t126 + (t146 / 0.2e1 - t298) * t130 + t311 * t97 + t99 * t10 + Ifges(4,5) * t95 + Ifges(4,6) * t96 + t313 * t47 + t54 * mrSges(4,1) - t55 * mrSges(4,2) - t45 * t56 - t44 * t57 - t41 * t58 - t38 * t59 - pkin(3) * t11 + (-t120 * t320 + t124 * t309 - t127 * t318 + t174) * g(3) + (-t21 * (-t127 * mrSges(6,1) - mrSges(6,2) * t237) - t30 * (mrSges(5,1) * t127 + mrSges(5,3) * t237) - t31 * (-mrSges(5,2) * t127 + mrSges(5,3) * t240) - t22 * (mrSges(6,2) * t240 + t127 * mrSges(6,3)) + (t134 / 0.2e1 - t133 / 0.2e1) * t88) * qJD(1) - t3 * t171 - t50 * t173 + (-t163 / 0.2e1 + t156 / 0.2e1) * t242 + t324 * t226 / 0.2e1 + ((-qJD(4) * t152 + t177) * m(6) + (-qJD(4) * t151 + t176) * m(5) + t320 * t124 * t265 + t256 * t219 - t257 * t220 + t325) * pkin(7) + (-t295 + t326) * mrSges(6,2) + (-t295 + t327) * mrSges(5,3) + (-m(6) * t228 + (-m(6) * t154 + t328) * t233) * g(1) + t323 * t225 + (t258 * t305 + t259 * t308 - t123 * t27 / 0.2e1 + t170 * t23 + t172 * t77) * t108 + t24 * t220 / 0.2e1 + Ifges(4,3) * qJDD(3) + (t258 * t303 + t259 * t307) * t89 - t101 * t241 + t155 * t283 + t162 * t284 - t175 * t266 + (-t21 * t41 - t22 * t38 + t313 * t23 + t3 * t99) * m(6) + (-pkin(3) * t50 - g(1) * t228 - t30 * t44 - t31 * t45 - t77 * t97) * m(5) + (t124 * t318 + t309 * t127 + t175) * t265; (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t22 - t23 * t46) * m(6) + t300 + (-t288 * t72 + t289 * t71) * g(1) + (t288 * t74 - t289 * t73) * g(2) + (-t317 * t88 - t331 * t89) * t333 + (-t332 * t88 + t24 - t260 + t82) * t334 + (-(-m(6) * t153 - t170) * t127 + t81) * g(3) - t77 * (mrSges(5,1) * t89 - mrSges(5,2) * t88) - t23 * (mrSges(6,1) * t89 + mrSges(6,3) * t88) + (-m(6) * t22 - t257 - t271) * t30 + (-m(6) * t21 - t256 + t270) * t31 + qJD(5) * t59 - t46 * t47 - pkin(4) * t15 + qJ(5) * t17 + t286 + (-Ifges(5,2) * t89 + t314 - t83) * t280 + t22 * t272 + t21 * t273 + t27 * t278 + (Ifges(6,3) * t89 - t269) * t281; -t108 * t59 + t89 * t47 + (-g(1) * t71 + g(2) * t73 - g(3) * t239 - t22 * t108 + t23 * t89 + t2) * m(6) + t15;];
tau = t6;
