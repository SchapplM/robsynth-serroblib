% Calculate vector of inverse dynamics joint torques for
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:10
% EndTime: 2019-12-31 17:18:26
% DurationCPUTime: 9.52s
% Computational Cost: add. (1830->381), mult. (4269->525), div. (0->0), fcn. (2470->6), ass. (0->181)
t312 = Ifges(4,4) + Ifges(5,4);
t313 = Ifges(4,1) + Ifges(5,1);
t291 = -Ifges(5,5) - Ifges(4,5);
t311 = Ifges(4,2) + Ifges(5,2);
t310 = Ifges(4,6) + Ifges(5,6);
t117 = sin(qJ(3));
t118 = sin(qJ(2));
t210 = qJD(1) * t118;
t192 = t117 * t210;
t120 = cos(qJ(3));
t207 = qJD(2) * t120;
t87 = -t192 + t207;
t319 = -t87 / 0.2e1;
t190 = t120 * t210;
t88 = qJD(2) * t117 + t190;
t252 = t88 / 0.2e1;
t121 = cos(qJ(2));
t209 = qJD(1) * t121;
t299 = t209 - qJD(3);
t318 = -t299 / 0.2e1;
t317 = mrSges(4,3) + mrSges(5,3);
t201 = qJD(1) * qJD(2);
t91 = qJDD(1) * t121 - t118 * t201;
t82 = t91 * pkin(5);
t316 = t121 * t82;
t315 = t312 * t120;
t314 = t312 * t117;
t92 = qJDD(1) * t118 + t121 * t201;
t35 = qJD(3) * t87 + qJDD(2) * t117 + t120 * t92;
t258 = t35 / 0.2e1;
t36 = -qJD(3) * t88 + qJDD(2) * t120 - t117 * t92;
t257 = t36 / 0.2e1;
t309 = -Ifges(5,3) - Ifges(4,3);
t110 = Ifges(3,4) * t209;
t302 = t312 * t87;
t287 = t291 * t299 + t313 * t88 + t302;
t308 = Ifges(3,1) * t210 + Ifges(3,5) * qJD(2) + t120 * t287 + t110;
t107 = pkin(3) * t120 + pkin(2);
t156 = -mrSges(5,1) * t120 + mrSges(5,2) * t117;
t158 = -mrSges(4,1) * t120 + mrSges(4,2) * t117;
t307 = m(4) * pkin(2) + m(5) * t107 - t156 - t158;
t277 = -t117 * t310 - t120 * t291;
t275 = -t117 * t311 + t315;
t273 = t120 * t313 - t314;
t116 = -qJ(4) - pkin(6);
t306 = -m(4) * pkin(6) + m(5) * t116 - t317;
t230 = Ifges(3,4) * t118;
t150 = Ifges(3,2) * t121 + t230;
t305 = t309 * t318 + t291 * t252 + t310 * t319 + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t150 / 0.2e1;
t301 = t312 * t88;
t288 = -t299 * t310 + t311 * t87 + t301;
t304 = -t288 / 0.2e1;
t84 = qJDD(3) - t91;
t303 = -t291 * t84 / 0.2e1 + t312 * t257 + t313 * t258;
t119 = sin(qJ(1));
t122 = cos(qJ(1));
t300 = g(1) * t122 + g(2) * t119;
t259 = m(5) * pkin(3);
t298 = t91 / 0.2e1;
t297 = t92 / 0.2e1;
t296 = t310 * t84 + t311 * t36 + t312 * t35;
t206 = qJD(2) * t121;
t194 = pkin(5) * t206;
t294 = -mrSges(5,1) - mrSges(4,1);
t293 = mrSges(2,2) - mrSges(3,3);
t292 = mrSges(4,2) + mrSges(5,2);
t170 = qJD(3) * t116;
t202 = qJD(4) * t120;
t215 = t118 * t120;
t217 = t117 * t121;
t161 = pkin(2) * t118 - pkin(6) * t121;
t89 = t161 * qJD(1);
t70 = t117 * t89;
t286 = t117 * t170 + t202 - t70 - (-pkin(5) * t215 - qJ(4) * t217) * qJD(1);
t213 = t120 * t121;
t139 = pkin(3) * t118 - qJ(4) * t213;
t48 = pkin(5) * t192 + t120 * t89;
t285 = -qJD(1) * t139 - qJD(4) * t117 + t120 * t170 - t48;
t284 = qJD(2) * mrSges(3,1) + mrSges(4,1) * t87 - mrSges(4,2) * t88 - mrSges(3,3) * t210;
t106 = pkin(5) * t213;
t162 = pkin(2) * t121 + pkin(6) * t118;
t94 = -pkin(1) - t162;
t56 = t117 * t94 + t106;
t112 = pkin(5) * t209;
t191 = t117 * t209;
t205 = qJD(3) * t117;
t283 = -t112 + (-t191 + t205) * pkin(3);
t282 = t259 + mrSges(5,1);
t281 = -t118 * t309 + t121 * t277;
t280 = t118 * t310 + t121 * t275;
t279 = -t118 * t291 + t121 * t273;
t276 = -t117 * t291 + t120 * t310;
t274 = t120 * t311 + t314;
t272 = t117 * t313 + t315;
t271 = -t291 * t35 - t309 * t84 + t310 * t36;
t83 = t92 * pkin(5);
t270 = t118 * t83 + t316;
t101 = qJD(2) * pkin(6) + t112;
t203 = qJD(3) * t120;
t219 = qJDD(1) * pkin(1);
t43 = -pkin(2) * t91 - pkin(6) * t92 - t219;
t62 = qJDD(2) * pkin(6) + t82;
t76 = t94 * qJD(1);
t7 = -t101 * t205 + t117 * t43 + t120 * t62 + t203 * t76;
t42 = t101 * t120 + t117 * t76;
t8 = -qJD(3) * t42 - t117 * t62 + t120 * t43;
t269 = -t117 * t8 + t120 * t7;
t22 = qJ(4) * t87 + t42;
t268 = -mrSges(4,3) * t42 - mrSges(5,3) * t22;
t41 = -t101 * t117 + t120 * t76;
t21 = -qJ(4) * t88 + t41;
t14 = -pkin(3) * t299 + t21;
t267 = -mrSges(4,3) * t41 - mrSges(5,3) * t14;
t266 = -m(5) - m(4) - m(3);
t265 = mrSges(4,1) + t282;
t160 = mrSges(3,1) * t121 - mrSges(3,2) * t118;
t262 = t118 * t317 + mrSges(2,1) + t160;
t2 = qJ(4) * t36 + qJD(4) * t87 + t7;
t261 = -t8 * mrSges(4,1) + t7 * mrSges(4,2) + t2 * mrSges(5,2);
t256 = t84 / 0.2e1;
t251 = pkin(3) * t88;
t245 = pkin(5) * t118;
t90 = t161 * qJD(2);
t232 = t117 * t90 + t203 * t94;
t231 = mrSges(5,2) * t120;
t229 = Ifges(3,4) * t121;
t208 = qJD(2) * t118;
t195 = pkin(5) * t208;
t220 = t117 * t195 + t120 * t90;
t218 = t117 * t118;
t216 = t117 * t122;
t214 = t119 * t117;
t212 = t121 * t122;
t204 = qJD(3) * t118;
t111 = pkin(5) * t210;
t189 = t117 * t206;
t10 = -mrSges(5,1) * t36 + mrSges(5,2) * t35;
t100 = -qJD(2) * pkin(2) + t111;
t63 = -qJDD(2) * pkin(2) + t83;
t159 = mrSges(3,1) * t118 + mrSges(3,2) * t121;
t157 = mrSges(4,1) * t117 + mrSges(4,2) * t120;
t155 = mrSges(5,1) * t117 + t231;
t145 = Ifges(3,5) * t121 - Ifges(3,6) * t118;
t140 = t107 * t121 - t116 * t118;
t138 = pkin(1) * t159;
t68 = -t117 * t212 + t119 * t120;
t66 = t120 * t122 + t121 * t214;
t135 = t118 * (Ifges(3,1) * t121 - t230);
t132 = -t117 * t204 + t120 * t206;
t131 = t118 * t203 + t189;
t99 = t116 * t120;
t97 = t116 * t117;
t96 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t209;
t93 = (pkin(3) * t117 + pkin(5)) * t118;
t86 = t120 * t94;
t77 = t157 * t118;
t69 = t120 * t212 + t214;
t67 = -t119 * t213 + t216;
t55 = -pkin(5) * t217 + t86;
t54 = pkin(3) * t131 + t194;
t53 = -mrSges(4,1) * t299 - mrSges(4,3) * t88;
t52 = -mrSges(5,1) * t299 - mrSges(5,3) * t88;
t51 = mrSges(4,2) * t299 + mrSges(4,3) * t87;
t50 = mrSges(5,2) * t299 + mrSges(5,3) * t87;
t49 = -pkin(5) * t190 + t70;
t47 = -pkin(3) * t87 + qJD(4) + t100;
t46 = -qJ(4) * t218 + t56;
t44 = -mrSges(5,1) * t87 + mrSges(5,2) * t88;
t40 = -qJ(4) * t215 + t86 + (-pkin(5) * t117 - pkin(3)) * t121;
t20 = -qJD(3) * t56 + t220;
t19 = (-t118 * t207 - t121 * t205) * pkin(5) + t232;
t18 = -mrSges(4,2) * t84 + mrSges(4,3) * t36;
t17 = -mrSges(5,2) * t84 + mrSges(5,3) * t36;
t16 = mrSges(4,1) * t84 - mrSges(4,3) * t35;
t15 = mrSges(5,1) * t84 - mrSges(5,3) * t35;
t13 = -pkin(3) * t36 + qJDD(4) + t63;
t12 = (-pkin(5) * qJD(2) - qJ(4) * qJD(3)) * t215 + (-qJD(4) * t118 + (-pkin(5) * qJD(3) - qJ(4) * qJD(2)) * t121) * t117 + t232;
t11 = -mrSges(4,1) * t36 + mrSges(4,2) * t35;
t9 = -t118 * t202 + t139 * qJD(2) + (-t106 + (qJ(4) * t118 - t94) * t117) * qJD(3) + t220;
t1 = pkin(3) * t84 - qJ(4) * t35 - qJD(4) * t88 + t8;
t3 = [t308 * t206 / 0.2e1 + (t41 * mrSges(4,1) + t14 * mrSges(5,1) - t42 * mrSges(4,2) - t22 * mrSges(5,2) - t305) * t208 + (t245 * t92 + t270 + t316) * mrSges(3,3) + t215 * t303 + t189 * t304 + t229 * t297 + t150 * t298 + t11 * t245 + (-t1 * mrSges(5,1) + Ifges(3,4) * t297 + Ifges(3,2) * t298 + t256 * t309 - t257 * t310 + t291 * t258 + t261) * t121 - (t117 * t287 + t120 * t288) * t204 / 0.2e1 + m(5) * (t1 * t40 + t12 * t22 + t13 * t93 + t14 * t9 + t2 * t46 + t47 * t54) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(5) * t270) - t271 * t121 / 0.2e1 + t100 * (mrSges(4,1) * t131 + mrSges(4,2) * t132) + t47 * (mrSges(5,1) * t131 + mrSges(5,2) * t132) - t284 * t194 + (t121 * (-Ifges(3,2) * t118 + t229) + t135) * t201 / 0.2e1 + (-t131 * t42 - t132 * t41 - t215 * t8 - t218 * t7) * mrSges(4,3) + t160 * t219 - pkin(1) * (-mrSges(3,1) * t91 + mrSges(3,2) * t92) + t93 * t10 + t63 * t77 + t54 * t44 + t55 * t16 + t56 * t18 + t46 * t17 + t12 * t50 + t19 * t51 + t9 * t52 + t20 * t53 + t40 * t15 + (qJD(2) * t279 - t204 * t272) * t252 + (qJD(2) * t280 - t204 * t274) * t87 / 0.2e1 - t138 * t201 - t96 * t195 + (qJD(2) * t281 - t204 * t276) * t318 + (-t1 * t215 - t131 * t22 - t132 * t14 - t2 * t218) * mrSges(5,3) + (-mrSges(3,1) * t245 + Ifges(3,5) * t118 + (-mrSges(3,2) * pkin(5) + Ifges(3,6)) * t121) * qJDD(2) + m(4) * (t100 * t194 + t19 * t42 + t20 * t41 + t55 * t8 + t56 * t7) - t296 * t218 / 0.2e1 + (m(4) * t63 * pkin(5) + Ifges(3,1) * t92 + Ifges(3,4) * t298 + t13 * t155 + t277 * t256 + t275 * t257 + t273 * t258) * t118 + (-t214 * t259 + t294 * t69 - t292 * t68 + t266 * (pkin(1) * t122 + pkin(5) * t119) + t293 * t119 + (-m(4) * t162 - m(5) * t140 - t262) * t122) * g(2) + (-t216 * t259 + t294 * t67 - t292 * t66 + (-m(5) * (-pkin(1) - t140) - m(4) * t94 + m(3) * pkin(1) + t262) * t119 + (pkin(5) * t266 + t293) * t122) * g(1) + qJD(2) ^ 2 * t145 / 0.2e1 + Ifges(2,3) * qJDD(1); -(-Ifges(3,2) * t210 + t110 + t308) * t209 / 0.2e1 + t300 * (t118 * t307 + t121 * t306 + t159) + (t118 * t306 - t121 * t307 - t160) * g(3) + t305 * t210 + t117 * t303 + (-t14 * (mrSges(5,1) * t118 - mrSges(5,3) * t213) - t41 * (mrSges(4,1) * t118 - mrSges(4,3) * t213) - t22 * (-mrSges(5,2) * t118 - mrSges(5,3) * t217) - t42 * (-mrSges(4,2) * t118 - mrSges(4,3) * t217) + (-t135 / 0.2e1 + t138) * qJD(1)) * qJD(1) + t269 * mrSges(4,3) + t299 * (-t100 * t157 - t47 * t155) + t284 * t112 + t285 * t52 + (t1 * t97 - t107 * t13 + t14 * t285 - t2 * t99 + t22 * t286 + t283 * t47) * m(5) + t286 * t50 + t288 * t191 / 0.2e1 + (t267 - t53 * pkin(6) + t287 / 0.2e1) * t203 - t107 * t10 + (t273 * t88 + t275 * t87 - t277 * t299) * qJD(3) / 0.2e1 - (t279 * t88 + t280 * t87 - t281 * t299) * qJD(1) / 0.2e1 + (-t1 * t117 + t120 * t2) * mrSges(5,3) + t96 * t111 + (t120 * t18 - t117 * t16 + m(4) * ((-t117 * t42 - t120 * t41) * qJD(3) + t269)) * pkin(6) + (-pkin(2) * t63 - t100 * t112 - t41 * t48 - t42 * t49) * m(4) + Ifges(3,5) * t92 + t97 * t15 - t99 * t17 + Ifges(3,6) * t91 - t82 * mrSges(3,2) - t83 * mrSges(3,1) - t49 * t51 - t48 * t53 - pkin(2) * t11 + t272 * t258 + t274 * t257 + t276 * t256 - t145 * t201 / 0.2e1 + (-pkin(6) * t51 + t268 + t304) * t205 + t13 * t156 + t63 * t158 + t283 * t44 + t296 * t120 / 0.2e1 + Ifges(3,3) * qJDD(2); t288 * t252 - t44 * t251 - m(5) * (t47 * t251 + (-t14 + t21) * t22) + (-mrSges(4,1) * t100 - mrSges(5,1) * t47 - t268) * t88 + (-t291 * t87 - t310 * t88) * t299 / 0.2e1 - (t313 * t87 - t301) * t88 / 0.2e1 + (-t311 * t88 + t287 + t302) * t319 + (t77 - (-t117 * t282 - t231) * t118) * g(3) + (-t265 * t68 + t292 * t69) * g(1) - t21 * t50 - t41 * t51 + t22 * t52 + t42 * t53 + pkin(3) * t15 + t282 * t1 + (t265 * t66 - t292 * t67) * g(2) + (-mrSges(4,2) * t100 - mrSges(5,2) * t47 - t267) * t87 - t261 + t271; -t87 * t50 + t88 * t52 + (g(3) * t121 - t118 * t300 + t14 * t88 - t22 * t87 + t13) * m(5) + t10;];
tau = t3;
