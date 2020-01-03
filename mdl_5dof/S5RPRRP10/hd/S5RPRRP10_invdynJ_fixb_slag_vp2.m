% Calculate vector of inverse dynamics joint torques for
% S5RPRRP10
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP10_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP10_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:55
% EndTime: 2019-12-31 18:51:20
% DurationCPUTime: 13.97s
% Computational Cost: add. (4366->459), mult. (10589->582), div. (0->0), fcn. (7545->10), ass. (0->203)
t337 = Ifges(5,4) + Ifges(6,4);
t150 = sin(qJ(4));
t153 = cos(qJ(4));
t146 = sin(pkin(8));
t147 = cos(pkin(8));
t151 = sin(qJ(3));
t154 = cos(qJ(3));
t124 = t146 * t151 - t154 * t147;
t121 = t124 * qJD(3);
t125 = t146 * t154 + t147 * t151;
t89 = -qJD(1) * t121 + qJDD(1) * t125;
t120 = t125 * qJD(1);
t98 = qJD(3) * t153 - t120 * t150;
t53 = qJD(4) * t98 + qJDD(3) * t150 + t153 * t89;
t277 = t53 / 0.2e1;
t99 = qJD(3) * t150 + t120 * t153;
t54 = -qJD(4) * t99 + qJDD(3) * t153 - t150 * t89;
t276 = t54 / 0.2e1;
t122 = t125 * qJD(3);
t90 = -qJD(1) * t122 - qJDD(1) * t124;
t83 = qJDD(4) - t90;
t275 = t83 / 0.2e1;
t338 = mrSges(5,3) + mrSges(6,3);
t304 = Ifges(5,1) + Ifges(6,1);
t303 = Ifges(5,5) + Ifges(6,5);
t302 = Ifges(5,2) + Ifges(6,2);
t301 = Ifges(5,6) + Ifges(6,6);
t318 = Ifges(6,3) + Ifges(5,3);
t217 = qJD(1) * qJD(2);
t135 = qJDD(1) * qJ(2) + t217;
t336 = t303 * t275 + t337 * t276 + t304 * t277;
t335 = t337 * t98;
t222 = t146 ^ 2 + t147 ^ 2;
t139 = pkin(4) * t153 + pkin(3);
t177 = -mrSges(6,1) * t153 + mrSges(6,2) * t150;
t179 = -mrSges(5,1) * t153 + mrSges(5,2) * t150;
t334 = -m(5) * pkin(3) - m(6) * t139 + t177 + t179;
t148 = -qJ(5) - pkin(7);
t333 = -m(5) * pkin(7) + m(6) * t148 - t338;
t332 = t337 * t153;
t331 = t337 * t150;
t330 = t337 * t99;
t242 = mrSges(4,3) * t120;
t298 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t98 - mrSges(5,2) * t99 - t242;
t247 = pkin(6) + qJ(2);
t130 = t247 * t146;
t126 = qJD(1) * t130;
t131 = t247 * t147;
t127 = qJD(1) * t131;
t91 = -t126 * t154 - t151 * t127;
t84 = -qJD(3) * pkin(3) - t91;
t329 = -m(5) * t84 + t298;
t138 = pkin(2) * t147 + pkin(1);
t129 = -qJD(1) * t138 + qJD(2);
t293 = Ifges(4,5) * qJD(3);
t328 = -t129 * mrSges(4,2) - t293 / 0.2e1;
t119 = t124 * qJD(1);
t66 = pkin(3) * t119 - pkin(7) * t120 + t129;
t92 = -t126 * t151 + t127 * t154;
t85 = qJD(3) * pkin(7) + t92;
t26 = -t150 * t85 + t153 * t66;
t17 = -qJ(5) * t99 + t26;
t312 = qJD(4) + t119;
t16 = pkin(4) * t312 + t17;
t27 = t150 * t66 + t153 * t85;
t18 = qJ(5) * t98 + t27;
t292 = Ifges(4,6) * qJD(3);
t327 = -t129 * mrSges(4,1) - t26 * mrSges(5,1) - t16 * mrSges(6,1) + t27 * mrSges(5,2) + t18 * mrSges(6,2) + t292 / 0.2e1;
t326 = m(6) * pkin(4);
t325 = t301 * t83 + t302 * t54 + t337 * t53;
t299 = t303 * t312 + t304 * t99 + t335;
t322 = t299 / 0.2e1;
t320 = -mrSges(6,1) - mrSges(5,1);
t305 = mrSges(5,2) + mrSges(6,2);
t317 = t301 * t98 + t303 * t99 + t318 * t312;
t300 = t301 * t312 + t302 * t98 + t330;
t107 = Ifges(4,4) * t119;
t316 = t119 * Ifges(4,2);
t289 = -t301 * t150 + t303 * t153;
t288 = -t302 * t150 + t332;
t287 = t304 * t153 - t331;
t314 = t301 * t54 + t303 * t53 + t318 * t83;
t152 = sin(qJ(1));
t155 = cos(qJ(1));
t313 = g(1) * t155 + g(2) * t152;
t311 = -m(4) - m(6) - m(5);
t196 = m(3) * qJ(2) + mrSges(3,3);
t310 = -t196 - mrSges(4,3) + mrSges(2,2);
t145 = pkin(8) + qJ(3);
t141 = sin(t145);
t142 = cos(t145);
t181 = mrSges(4,1) * t142 - mrSges(4,2) * t141;
t182 = -mrSges(3,1) * t147 + mrSges(3,2) * t146;
t308 = m(3) * pkin(1) + t338 * t141 + mrSges(2,1) + t181 - t182;
t218 = qJD(4) * t153;
t219 = qJD(4) * t150;
t128 = -qJDD(1) * t138 + qJDD(2);
t33 = -pkin(3) * t90 - pkin(7) * t89 + t128;
t188 = pkin(6) * qJDD(1) + t135;
t105 = t188 * t146;
t106 = t188 * t147;
t220 = qJD(3) * t154;
t221 = qJD(3) * t151;
t36 = -t151 * t105 + t154 * t106 - t126 * t220 - t127 * t221;
t34 = qJDD(3) * pkin(7) + t36;
t3 = t150 * t33 + t153 * t34 + t66 * t218 - t219 * t85;
t2 = qJ(5) * t54 + qJD(5) * t98 + t3;
t4 = -qJD(4) * t27 - t150 * t34 + t153 * t33;
t307 = t4 * mrSges(5,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2);
t189 = qJD(4) * t148;
t232 = t119 * t153;
t86 = pkin(3) * t120 + pkin(7) * t119;
t38 = -t150 * t91 + t153 * t86;
t297 = -pkin(4) * t120 - qJ(5) * t232 - qJD(5) * t150 + t153 * t189 - t38;
t233 = t119 * t150;
t296 = -t92 + (t219 + t233) * pkin(4);
t39 = t150 * t86 + t153 * t91;
t295 = -qJ(5) * t233 + qJD(5) * t153 + t150 * t189 - t39;
t294 = mrSges(6,1) + t326;
t291 = -t154 * t130 - t131 * t151;
t244 = mrSges(6,2) * t153;
t176 = mrSges(6,1) * t150 + t244;
t178 = mrSges(5,1) * t150 + mrSges(5,2) * t153;
t55 = -pkin(4) * t98 + qJD(5) + t84;
t290 = t55 * t176 + t84 * t178;
t286 = -t27 * mrSges(5,3) - t18 * mrSges(6,3);
t285 = -t26 * mrSges(5,3) - t16 * mrSges(6,3);
t284 = -t150 * t4 + t153 * t3;
t280 = mrSges(5,1) + t294;
t274 = -t98 / 0.2e1;
t273 = t98 / 0.2e1;
t272 = -t99 / 0.2e1;
t271 = t99 / 0.2e1;
t270 = pkin(4) * t99;
t269 = -t312 / 0.2e1;
t268 = t312 / 0.2e1;
t267 = t119 / 0.2e1;
t265 = t120 / 0.2e1;
t88 = pkin(3) * t124 - pkin(7) * t125 - t138;
t95 = -t130 * t151 + t131 * t154;
t93 = t153 * t95;
t52 = t150 * t88 + t93;
t243 = mrSges(4,3) * t119;
t241 = Ifges(4,4) * t120;
t231 = t121 * t153;
t230 = t125 * t150;
t229 = t125 * t153;
t227 = t150 * t121;
t226 = t150 * t152;
t225 = t150 * t155;
t224 = t152 * t153;
t223 = t153 * t155;
t215 = qJDD(1) * t146;
t214 = qJDD(1) * t147;
t67 = -t124 * qJD(2) + qJD(3) * t291;
t87 = pkin(3) * t122 + pkin(7) * t121;
t211 = t150 * t87 + t153 * t67 + t88 * t218;
t205 = t125 * t218;
t194 = -t90 * mrSges(4,1) + t89 * mrSges(4,2);
t14 = -t54 * mrSges(6,1) + t53 * mrSges(6,2);
t192 = -t219 / 0.2e1;
t190 = -t150 * t67 + t153 * t87;
t51 = -t150 * t95 + t153 * t88;
t184 = pkin(3) * t142 + pkin(7) * t141;
t183 = -mrSges(3,1) * t214 + mrSges(3,2) * t215;
t167 = t139 * t142 - t141 * t148;
t166 = qJ(5) * t121 - qJD(5) * t125;
t37 = -t105 * t154 - t151 * t106 + t126 * t221 - t127 * t220;
t114 = -t142 * t225 + t224;
t112 = t142 * t226 + t223;
t163 = t205 - t227;
t162 = t125 * t219 + t231;
t35 = -qJDD(3) * pkin(3) - t37;
t68 = qJD(2) * t125 + qJD(3) * t95;
t140 = -qJDD(1) * pkin(1) + qJDD(2);
t133 = t148 * t153;
t132 = t148 * t150;
t115 = t142 * t223 + t226;
t113 = -t142 * t224 + t225;
t100 = -qJD(3) * mrSges(4,2) - t243;
t72 = t120 * Ifges(4,1) - t107 + t293;
t71 = t241 + t292 - t316;
t69 = pkin(4) * t230 - t291;
t65 = mrSges(5,1) * t312 - mrSges(5,3) * t99;
t64 = mrSges(6,1) * t312 - mrSges(6,3) * t99;
t63 = -mrSges(5,2) * t312 + mrSges(5,3) * t98;
t62 = -mrSges(6,2) * t312 + mrSges(6,3) * t98;
t56 = -mrSges(6,1) * t98 + mrSges(6,2) * t99;
t32 = pkin(4) * t163 + t68;
t28 = -qJ(5) * t230 + t52;
t24 = pkin(4) * t124 - qJ(5) * t229 + t51;
t23 = -mrSges(5,2) * t83 + mrSges(5,3) * t54;
t22 = -mrSges(6,2) * t83 + mrSges(6,3) * t54;
t21 = mrSges(5,1) * t83 - mrSges(5,3) * t53;
t20 = mrSges(6,1) * t83 - mrSges(6,3) * t53;
t15 = -mrSges(5,1) * t54 + mrSges(5,2) * t53;
t13 = -pkin(4) * t54 + qJDD(5) + t35;
t12 = -qJD(4) * t52 + t190;
t11 = -t219 * t95 + t211;
t6 = -qJ(5) * t205 + (-qJD(4) * t95 + t166) * t150 + t211;
t5 = pkin(4) * t122 + t166 * t153 + (-t93 + (qJ(5) * t125 - t88) * t150) * qJD(4) + t190;
t1 = pkin(4) * t83 - qJ(5) * t53 - qJD(5) * t99 + t4;
t7 = [(-t162 * t303 - t163 * t301) * t268 - t299 * t231 / 0.2e1 + (t162 * t26 - t163 * t27 - t229 * t4 - t230 * t3) * mrSges(5,3) + (t128 * mrSges(4,1) + t1 * mrSges(6,1) - t36 * mrSges(4,3) - Ifges(4,4) * t89 - Ifges(4,2) * t90 - Ifges(4,6) * qJDD(3) + t318 * t275 + t301 * t276 + t303 * t277 + t307 + t314 / 0.2e1) * t124 + (t128 * mrSges(4,2) - t37 * mrSges(4,3) + Ifges(4,1) * t89 + Ifges(4,4) * t90 + Ifges(4,5) * qJDD(3) + t13 * t176 + t178 * t35 + t192 * t299 + t275 * t289 + t276 * t288 + t277 * t287) * t125 + t140 * t182 - pkin(1) * t183 - t138 * t194 - t325 * t230 / 0.2e1 + (-t226 * t326 + t311 * (t155 * t138 + t152 * t247) + t320 * t115 - t305 * t114 + t310 * t152 + (-m(5) * t184 - m(6) * t167 - t308) * t155) * g(2) + (t320 * t113 - t305 * t112 + (-t150 * t326 + t247 * t311 + t310) * t155 + (m(4) * t138 - m(6) * (-t138 - t167) - m(5) * (-t138 - t184) + t308) * t152) * g(1) + t84 * (mrSges(5,1) * t163 - mrSges(5,2) * t162) + t55 * (mrSges(6,1) * t163 - mrSges(6,2) * t162) + (-t205 / 0.2e1 + t227 / 0.2e1) * t300 + (-t162 * t337 - t163 * t302) * t273 + (-t162 * t304 - t163 * t337) * t271 + 0.2e1 * t222 * t135 * mrSges(3,3) + t229 * t336 + (-t1 * t229 + t16 * t162 - t163 * t18 - t2 * t230) * mrSges(6,3) + (Ifges(3,4) * t146 + Ifges(3,2) * t147) * t214 + (Ifges(3,1) * t146 + Ifges(3,4) * t147) * t215 + m(5) * (t11 * t27 + t12 * t26 - t291 * t35 + t3 * t52 + t4 * t51) + m(4) * (-t128 * t138 + t291 * t37 + t36 * t95 + t67 * t92) + (mrSges(4,1) * t291 - mrSges(4,2) * t95) * qJDD(3) + (t121 * t91 - t291 * t89 + t90 * t95) * mrSges(4,3) - t291 * t15 + (-m(4) * t91 - t329) * t68 + (t317 / 0.2e1 - Ifges(4,4) * t265 + t316 / 0.2e1 - t71 / 0.2e1 + t301 * t273 + t303 * t271 + t318 * t268 - t92 * mrSges(4,3) - t327) * t122 + (-Ifges(4,1) * t265 + t107 / 0.2e1 - t72 / 0.2e1 + t328) * t121 + t24 * t20 + t28 * t22 + t51 * t21 + t52 * t23 + t32 * t56 + t6 * t62 + t11 * t63 + t5 * t64 + t12 * t65 + t69 * t14 + t67 * t100 + m(3) * (-pkin(1) * t140 + (t135 + t217) * qJ(2) * t222) + m(6) * (t1 * t24 + t13 * t69 + t16 * t5 + t18 * t6 + t2 * t28 + t32 * t55) + Ifges(2,3) * qJDD(1); t119 * t100 + (-t56 + t298) * t120 + (t20 + t21 + t312 * (t62 + t63)) * t153 + (t22 + t23 - t312 * (t64 + t65)) * t150 + m(3) * t140 + t183 + t194 + (-g(1) * t152 + g(2) * t155) * (m(3) - t311) - t196 * t222 * qJD(1) ^ 2 + (t1 * t153 - t120 * t55 + t150 * t2 + t312 * (-t150 * t16 + t153 * t18)) * m(6) + (-t120 * t84 + t150 * t3 + t153 * t4 + t312 * (-t150 * t26 + t153 * t27)) * m(5) + (t119 * t92 + t120 * t91 + t128) * m(4); (-t65 * pkin(7) + t285 + t322) * t218 - (-Ifges(4,1) * t119 - t241 + t317) * t120 / 0.2e1 + (-t233 / 0.2e1 + t192) * t300 + t232 * t322 + (-t1 * t150 + t153 * t2 - t16 * t232 - t18 * t233) * mrSges(6,3) + (-t243 - t100) * t91 + (t72 - t107) * t267 + (-t232 * t26 - t233 * t27 + t284) * mrSges(5,3) + t286 * t219 + t290 * qJD(4) + t13 * t177 + t35 * t179 - g(3) * t181 + (-t63 * t219 - t150 * t21 + t153 * t23 + m(5) * ((-t150 * t27 - t153 * t26) * qJD(4) + t284)) * pkin(7) + t325 * t153 / 0.2e1 + t150 * t336 + t71 * t265 + (t287 * t99 + t288 * t98 + t289 * t312) * qJD(4) / 0.2e1 + (t242 + t329) * t92 + (t153 * t302 + t331) * t276 + (t150 * t304 + t332) * t277 + (t313 * (mrSges(4,1) - t334) + t333 * g(3)) * t141 + (t313 * (mrSges(4,2) + t333) + t334 * g(3)) * t142 + (-Ifges(4,2) * t267 + t318 * t269 + t303 * t272 + t301 * t274 + t327) * t120 - (t289 * t269 + t287 * t272 + t288 * t274 - t290 + t328) * t119 + (-pkin(3) * t35 - t26 * t38 - t27 * t39) * m(5) - pkin(3) * t15 + t295 * t62 + t296 * t56 + t297 * t64 + (t1 * t132 - t13 * t139 - t133 * t2 + t16 * t297 + t18 * t295 + t296 * t55) * m(6) - t36 * mrSges(4,2) + t37 * mrSges(4,1) + (t150 * t303 + t153 * t301) * t275 - t39 * t63 - t38 * t65 + Ifges(4,5) * t89 + Ifges(4,6) * t90 + Ifges(4,3) * qJDD(3) + t132 * t20 - t133 * t22 - t139 * t14; t307 + (t150 * t294 + t178 + t244) * g(3) * t141 + (-t84 * mrSges(5,2) - t55 * mrSges(6,2) - t285) * t98 + (-t84 * mrSges(5,1) - t55 * mrSges(6,1) - t286) * t99 + t294 * t1 + pkin(4) * t20 + t300 * t271 + (-t302 * t99 + t299 + t335) * t274 + (-t301 * t99 + t303 * t98) * t269 + (t304 * t98 - t330) * t272 + (-t280 * t114 + t115 * t305) * g(1) + (t280 * t112 - t113 * t305) * g(2) - t17 * t62 - t26 * t63 + t18 * t64 + t27 * t65 - t56 * t270 - m(6) * (t55 * t270 + (-t16 + t17) * t18) + t314; -t98 * t62 + t99 * t64 + (g(3) * t142 - t141 * t313 + t16 * t99 - t18 * t98 + t13) * m(6) + t14;];
tau = t7;
