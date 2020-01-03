% Calculate vector of inverse dynamics joint torques for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP7_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:10
% EndTime: 2019-12-31 17:20:24
% DurationCPUTime: 8.95s
% Computational Cost: add. (1845->402), mult. (4261->548), div. (0->0), fcn. (2451->6), ass. (0->191)
t303 = Ifges(4,1) + Ifges(5,1);
t302 = Ifges(5,4) + Ifges(4,5);
t301 = Ifges(4,6) - Ifges(5,6);
t116 = sin(qJ(3));
t117 = sin(qJ(2));
t202 = qJD(1) * t117;
t184 = t116 * t202;
t119 = cos(qJ(3));
t194 = t119 * qJD(2);
t82 = t184 - t194;
t250 = t82 / 0.2e1;
t183 = t119 * t202;
t83 = qJD(2) * t116 + t183;
t306 = -t83 / 0.2e1;
t120 = cos(qJ(2));
t201 = qJD(1) * t120;
t96 = qJD(3) - t201;
t305 = t96 / 0.2e1;
t193 = qJD(1) * qJD(2);
t86 = qJDD(1) * t120 - t117 * t193;
t78 = t86 * pkin(5);
t304 = -mrSges(4,3) - mrSges(5,2);
t300 = -Ifges(4,3) - Ifges(5,2);
t276 = -t301 * t116 + t302 * t119;
t224 = Ifges(5,5) * t116;
t226 = Ifges(4,4) * t116;
t273 = t303 * t119 + t224 - t226;
t158 = t119 * mrSges(5,1) + t116 * mrSges(5,3);
t160 = mrSges(4,1) * t119 - mrSges(4,2) * t116;
t299 = -t158 - t160;
t106 = Ifges(3,4) * t201;
t76 = Ifges(5,5) * t83;
t24 = Ifges(5,6) * t96 + Ifges(5,3) * t82 + t76;
t242 = Ifges(5,5) * t82;
t77 = Ifges(4,4) * t82;
t287 = t302 * t96 + t303 * t83 + t242 - t77;
t298 = Ifges(3,1) * t202 + Ifges(3,5) * qJD(2) + t116 * t24 + t119 * t287 + t106;
t228 = Ifges(3,4) * t117;
t282 = t120 * Ifges(3,2);
t152 = t228 + t282;
t297 = Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t152 / 0.2e1 + t300 * t305 + t301 * t250 + t302 * t306;
t214 = qJD(3) * t82;
t87 = qJDD(1) * t117 + t120 * t193;
t35 = qJDD(2) * t116 + t119 * t87 - t214;
t255 = t35 / 0.2e1;
t36 = qJD(3) * t83 - t119 * qJDD(2) + t116 * t87;
t254 = -t36 / 0.2e1;
t80 = qJDD(3) - t86;
t252 = t80 / 0.2e1;
t296 = t86 / 0.2e1;
t295 = t87 / 0.2e1;
t294 = -m(4) - m(5);
t293 = t302 * t80 + (-Ifges(4,4) + Ifges(5,5)) * t36 + t303 * t35;
t199 = qJD(2) * t120;
t292 = pkin(5) * t199;
t291 = mrSges(2,2) - mrSges(3,3);
t204 = t120 * pkin(2) + t117 * pkin(6);
t289 = -pkin(1) - t204;
t50 = -mrSges(4,2) * t96 - mrSges(4,3) * t82;
t53 = -mrSges(5,2) * t82 + mrSges(5,3) * t96;
t286 = -t50 - t53;
t51 = mrSges(4,1) * t96 - mrSges(4,3) * t83;
t52 = -mrSges(5,1) * t96 + mrSges(5,2) * t83;
t285 = t51 - t52;
t108 = pkin(5) * t201;
t141 = pkin(3) * t116 - qJ(4) * t119;
t284 = -qJD(4) * t116 + t141 * t96 - t108;
t283 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t82 + mrSges(4,2) * t83 + mrSges(3,3) * t202;
t162 = mrSges(3,1) * t120 - mrSges(3,2) * t117;
t281 = -mrSges(2,1) - t162;
t279 = -t117 * t300 + t120 * t276;
t278 = t117 * t302 + t120 * t273;
t157 = t116 * mrSges(5,1) - t119 * mrSges(5,3);
t159 = mrSges(4,1) * t116 + mrSges(4,2) * t119;
t107 = pkin(5) * t202;
t93 = -qJD(2) * pkin(2) + t107;
t23 = pkin(3) * t82 - qJ(4) * t83 + t93;
t277 = t23 * t157 + t93 * t159;
t275 = t116 * t302 + t119 * t301;
t223 = Ifges(5,5) * t119;
t225 = Ifges(4,4) * t119;
t274 = t116 * t303 - t223 + t225;
t272 = -t300 * t80 - t301 * t36 + t302 * t35;
t79 = t87 * pkin(5);
t271 = t117 * t79 + t120 * t78;
t270 = t304 * t117;
t196 = qJD(3) * t119;
t198 = qJD(3) * t116;
t213 = qJDD(1) * pkin(1);
t40 = -pkin(2) * t86 - pkin(6) * t87 - t213;
t61 = qJDD(2) * pkin(6) + t78;
t73 = t289 * qJD(1);
t94 = qJD(2) * pkin(6) + t108;
t8 = t116 * t40 + t119 * t61 + t73 * t196 - t198 * t94;
t39 = t116 * t73 + t119 * t94;
t9 = -qJD(3) * t39 - t116 * t61 + t119 * t40;
t269 = -t9 * t116 + t8 * t119;
t22 = qJ(4) * t96 + t39;
t268 = -t22 * mrSges(5,2) - t39 * mrSges(4,3);
t38 = -t116 * t94 + t119 * t73;
t265 = -t38 + qJD(4);
t21 = -pkin(3) * t96 + t265;
t267 = t21 * mrSges(5,2) - t38 * mrSges(4,3);
t1 = qJ(4) * t80 + qJD(4) * t96 + t8;
t3 = -pkin(3) * t80 + qJDD(4) - t9;
t266 = t1 * t119 + t3 * t116;
t264 = -t35 * Ifges(5,5) / 0.2e1 - t80 * Ifges(5,6) / 0.2e1 + Ifges(4,4) * t255 + Ifges(4,6) * t252 + (Ifges(5,3) + Ifges(4,2)) * t254;
t263 = pkin(6) * t294;
t161 = mrSges(3,1) * t117 + mrSges(3,2) * t120;
t142 = pkin(3) * t119 + qJ(4) * t116;
t88 = -pkin(2) - t142;
t262 = t161 + t304 * t120 + (m(4) * pkin(2) - m(5) * t88 - t299) * t117;
t195 = qJD(3) * t120;
t200 = qJD(2) * t117;
t163 = pkin(2) * t117 - pkin(6) * t120;
t85 = t163 * qJD(2);
t20 = pkin(5) * (t116 * t200 - t119 * t195) + t119 * t85 - t198 * t289;
t260 = m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1);
t259 = m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3);
t258 = -t9 * mrSges(4,1) + t3 * mrSges(5,1) + t8 * mrSges(4,2) - t1 * mrSges(5,3);
t243 = Ifges(4,4) * t83;
t27 = -Ifges(4,2) * t82 + Ifges(4,6) * t96 + t243;
t256 = -t27 / 0.2e1;
t253 = t36 / 0.2e1;
t251 = -t82 / 0.2e1;
t248 = t83 / 0.2e1;
t245 = t116 / 0.2e1;
t241 = pkin(5) * t117;
t207 = t119 * t120;
t56 = pkin(5) * t207 + t116 * t289;
t227 = Ifges(3,4) * t120;
t84 = t163 * qJD(1);
t218 = t119 * t84;
t217 = t119 * t289;
t212 = t116 * t117;
t211 = t116 * t120;
t210 = t117 * t119;
t121 = cos(qJ(1));
t209 = t117 * t121;
t118 = sin(qJ(1));
t208 = t118 * t120;
t206 = t120 * t121;
t205 = t121 * t116;
t203 = t121 * pkin(1) + t118 * pkin(5);
t197 = qJD(3) * t117;
t185 = pkin(5) * t116 + pkin(3);
t182 = t116 * t199;
t14 = -t80 * mrSges(5,1) + t35 * mrSges(5,2);
t169 = t196 / 0.2e1;
t168 = t193 / 0.2e1;
t62 = -qJDD(2) * pkin(2) + t79;
t151 = -Ifges(4,2) * t116 + t225;
t150 = Ifges(4,2) * t119 + t226;
t147 = Ifges(3,5) * t120 - Ifges(3,6) * t117;
t144 = Ifges(5,3) * t116 + t223;
t143 = -Ifges(5,3) * t119 + t224;
t138 = pkin(5) + t141;
t136 = pkin(1) * t161;
t133 = t117 * (Ifges(3,1) * t120 - t228);
t132 = -t116 * t197 + t120 * t194;
t131 = t117 * t196 + t182;
t126 = Ifges(4,6) * t117 + t120 * t151;
t125 = Ifges(5,6) * t117 + t120 * t144;
t19 = t116 * t85 + t289 * t196 + (-t116 * t195 - t117 * t194) * pkin(5);
t91 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t201;
t74 = t159 * t117;
t69 = t116 * t84;
t68 = t118 * t116 + t119 * t206;
t67 = -t118 * t119 + t120 * t205;
t66 = t118 * t207 - t205;
t65 = t116 * t208 + t119 * t121;
t59 = t138 * t117;
t55 = -pkin(5) * t211 + t217;
t49 = -pkin(5) * t183 + t69;
t48 = pkin(5) * t184 + t218;
t47 = t120 * t185 - t217;
t46 = -qJ(4) * t120 + t56;
t44 = mrSges(5,1) * t82 - mrSges(5,3) * t83;
t43 = pkin(3) * t83 + qJ(4) * t82;
t42 = -t185 * t202 - t218;
t41 = t69 + (-pkin(5) * t119 + qJ(4)) * t202;
t18 = (qJD(3) * t142 - qJD(4) * t119) * t117 + t138 * t199;
t17 = -pkin(3) * t200 - t20;
t16 = -mrSges(5,2) * t36 + mrSges(5,3) * t80;
t15 = -mrSges(4,2) * t80 - mrSges(4,3) * t36;
t13 = mrSges(4,1) * t80 - mrSges(4,3) * t35;
t12 = qJ(4) * t200 - qJD(4) * t120 + t19;
t11 = mrSges(4,1) * t36 + mrSges(4,2) * t35;
t10 = mrSges(5,1) * t36 - mrSges(5,3) * t35;
t2 = pkin(3) * t36 - qJ(4) * t35 - qJD(4) * t83 + t62;
t4 = [((-mrSges(3,2) * pkin(5) + Ifges(3,6)) * qJDD(2) + t78 * mrSges(3,3) - t272 / 0.2e1 + Ifges(3,4) * t295 + Ifges(3,2) * t296 - Ifges(4,6) * t254 - Ifges(5,6) * t253 + t227 * t168 - t302 * t255 + t300 * t252 + t258) * t120 + m(5) * (t1 * t46 + t12 * t22 + t17 * t21 + t18 * t23 + t2 * t59 + t3 * t47) + (t241 * t87 + t271) * mrSges(3,3) + t93 * (mrSges(4,1) * t131 + mrSges(4,2) * t132) + t23 * (mrSges(5,1) * t131 - mrSges(5,3) * t132) - t136 * t193 + (-m(3) * t203 + t294 * (pkin(2) * t206 + pkin(6) * t209 + t203) - t260 * t68 - t259 * t67 + t304 * t209 + t281 * t121 + t291 * t118) * g(2) + (-mrSges(3,1) * t241 + Ifges(3,5) * t117) * qJDD(2) + t182 * t256 + (t38 * mrSges(4,1) - t21 * mrSges(5,1) - t39 * mrSges(4,2) + t22 * mrSges(5,3) - pkin(5) * t91 - t297) * t200 - (t116 * t287 + t119 * t27) * t197 / 0.2e1 + t11 * t241 + (qJD(2) * t125 - t143 * t197) * t250 + (qJD(2) * t126 - t150 * t197) * t251 - pkin(1) * (-mrSges(3,1) * t86 + mrSges(3,2) * t87) + t62 * t74 + t19 * t50 + t20 * t51 + t17 * t52 + t12 * t53 + t55 * t13 + t56 * t15 + t59 * t10 + t18 * t44 + t46 * t16 + t47 * t14 + m(4) * (t19 * t39 + t20 * t38 + t292 * t93 + t55 * t9 + t56 * t8) + (-t1 * mrSges(5,2) - t8 * mrSges(4,3) - t264) * t212 + t133 * t168 + (-t131 * t39 - t132 * t38 - t210 * t9) * mrSges(4,3) + qJD(2) ^ 2 * t147 / 0.2e1 + (-t131 * t22 + t132 * t21 + t210 * t3) * mrSges(5,2) + m(3) * (qJDD(1) * pkin(1) ^ 2 + t271 * pkin(5)) + (qJD(2) * t278 - t197 * t274) * t248 + t162 * t213 + (t260 * t66 + t259 * t65 + (m(3) * pkin(1) + t289 * t294 - t270 - t281) * t118 + (t291 + (-m(3) + t294) * pkin(5)) * t121) * g(1) + (qJD(2) * t279 - t197 * t275) * t305 + t298 * t199 / 0.2e1 + t283 * t292 + t227 * t295 + t152 * t296 + t293 * t210 / 0.2e1 + (m(4) * pkin(5) * t62 + Ifges(3,1) * t87 + Ifges(3,4) * t296 + t144 * t253 + t151 * t254 + t2 * t157 - t168 * t282 + t24 * t169 + t276 * t252 + t273 * t255) * t117 + Ifges(2,3) * qJDD(1); -t283 * t108 + t284 * t44 - t147 * t193 / 0.2e1 + (-t151 / 0.2e1 + t144 / 0.2e1) * t214 + t150 * t254 + t297 * t202 + t143 * t253 + Ifges(3,6) * t86 + Ifges(3,5) * t87 + t88 * t10 - t78 * mrSges(3,2) - t79 * mrSges(3,1) - t49 * t50 - t48 * t51 - t42 * t52 - t41 * t53 + t266 * mrSges(5,2) + t267 * t196 + (t256 + t24 / 0.2e1 + t268) * t198 + (t118 * t262 + t208 * t263) * g(2) + (t121 * t262 + t206 * t263) * g(1) + t264 * t119 + (t2 * t88 - t21 * t42 - t22 * t41 + t284 * t23) * m(5) + (-pkin(2) * t62 - t108 * t93 - t38 * t48 - t39 * t49) * m(4) + (-t285 * t196 + t286 * t198 + ((-t22 * t116 + t21 * t119) * qJD(3) + t266) * m(5) + ((-t39 * t116 - t38 * t119) * qJD(3) + t269) * m(4) + (t16 + t15) * t119 + (t14 - t13) * t116) * pkin(6) + t91 * t107 + t269 * mrSges(4,3) + t274 * t255 + t275 * t252 + (t245 * t27 - t277) * t201 + t277 * qJD(3) + t287 * t169 - t2 * t158 + (t273 * t83 + t276 * t96) * qJD(3) / 0.2e1 - (t278 * t83 + t279 * t96) * qJD(1) / 0.2e1 - pkin(2) * t11 - t62 * t160 - (-Ifges(3,2) * t202 + t106 + t298) * t201 / 0.2e1 + (-t162 + t294 * t204 + (-m(5) * t142 + t299) * t120 + t270) * g(3) + (-t38 * (mrSges(4,1) * t117 - mrSges(4,3) * t207) - t21 * (-mrSges(5,1) * t117 + mrSges(5,2) * t207) - t39 * (-mrSges(4,2) * t117 - mrSges(4,3) * t211) - t22 * (-mrSges(5,2) * t211 + mrSges(5,3) * t117) + (t126 / 0.2e1 - t125 / 0.2e1) * t82 + (-t133 / 0.2e1 + t136) * qJD(1)) * qJD(1) + t293 * t245 + Ifges(3,3) * qJDD(2); (-t259 * t66 + t260 * t65) * g(2) + (-mrSges(4,1) * t93 - mrSges(5,1) * t23 + Ifges(5,3) * t251 - t268) * t83 + (mrSges(4,2) * t93 - mrSges(5,3) * t23 + t267) * t82 + t285 * t39 - (-t301 * t83 - t302 * t82) * t96 / 0.2e1 + (-t303 * t82 + t24 - t243 + t76) * t306 + (t74 - (-m(5) * t141 - t157) * t117) * g(3) + t286 * t38 - t242 * t251 + t27 * t248 + (-Ifges(4,2) * t83 + t287 - t77) * t250 + qJD(4) * t53 - t43 * t44 + qJ(4) * t16 - pkin(3) * t14 + (-pkin(3) * t3 + qJ(4) * t1 - t21 * t39 + t22 * t265 - t23 * t43) * m(5) - t258 + (-t259 * t68 + t260 * t67) * g(1) + t272; t83 * t44 - t96 * t53 + (-g(1) * t67 - g(2) * t65 - g(3) * t212 - t22 * t96 + t23 * t83 + t3) * m(5) + t14;];
tau = t4;
