% Calculate vector of inverse dynamics joint torques for
% S5RRPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:01
% EndTime: 2019-12-31 19:54:21
% DurationCPUTime: 10.34s
% Computational Cost: add. (4609->472), mult. (10779->602), div. (0->0), fcn. (7624->12), ass. (0->210)
t322 = -mrSges(5,1) - mrSges(6,1);
t311 = Ifges(5,1) + Ifges(6,1);
t309 = Ifges(6,4) + Ifges(5,5);
t189 = sin(qJ(2));
t191 = cos(qJ(2));
t155 = -mrSges(3,1) * t191 + mrSges(3,2) * t189;
t184 = qJ(2) + pkin(8);
t177 = sin(t184);
t178 = cos(t184);
t179 = qJ(4) + t184;
t167 = sin(t179);
t168 = cos(t179);
t300 = t322 * t168 + (mrSges(5,2) - mrSges(6,3)) * t167;
t321 = -mrSges(4,1) * t178 + mrSges(4,2) * t177 + t155 + t300;
t183 = qJD(2) + qJD(4);
t188 = sin(qJ(4));
t185 = sin(pkin(8));
t186 = cos(pkin(8));
t135 = -t185 * t189 + t186 * t191;
t121 = t135 * qJD(1);
t268 = pkin(7) * t121;
t187 = -qJ(3) - pkin(6);
t154 = t187 * t189;
t141 = qJD(1) * t154;
t132 = qJD(2) * pkin(2) + t141;
t156 = t187 * t191;
t142 = qJD(1) * t156;
t243 = t186 * t142;
t83 = t185 * t132 - t243;
t55 = t83 + t268;
t250 = t188 * t55;
t277 = cos(qJ(4));
t235 = t191 * qJD(1);
t236 = t189 * qJD(1);
t122 = -t185 * t235 - t186 * t236;
t267 = pkin(7) * t122;
t125 = t185 * t142;
t82 = t186 * t132 + t125;
t53 = qJD(2) * pkin(3) + t267 + t82;
t18 = t277 * t53 - t250;
t304 = qJD(5) - t18;
t11 = -t183 * pkin(4) + t304;
t202 = t188 * t121 - t122 * t277;
t75 = t277 * t121 + t188 * t122;
t273 = Ifges(6,5) * t75;
t68 = Ifges(5,4) * t75;
t307 = t309 * t183 + t311 * t202 - t273 + t68;
t180 = t191 * pkin(2);
t170 = t180 + pkin(1);
t147 = -qJD(1) * t170 + qJD(3);
t93 = -pkin(3) * t121 + t147;
t31 = -pkin(4) * t75 - qJ(5) * t202 + t93;
t320 = -t93 * mrSges(5,2) - t11 * mrSges(6,2) + t18 * mrSges(5,3) + t31 * mrSges(6,3) - t307 / 0.2e1;
t310 = -Ifges(5,4) + Ifges(6,5);
t308 = -Ifges(5,6) + Ifges(6,6);
t271 = pkin(2) * t186;
t169 = pkin(3) + t271;
t272 = pkin(2) * t185;
t117 = t188 * t169 + t277 * t272;
t89 = -t141 * t185 + t243;
t206 = t89 - t268;
t90 = t186 * t141 + t125;
t57 = t90 + t267;
t305 = t117 * qJD(4) - t188 * t57 + t206 * t277;
t190 = sin(qJ(1));
t192 = cos(qJ(1));
t318 = g(1) * t192 + g(2) * t190;
t317 = t168 * (-m(6) * qJ(5) - mrSges(6,3));
t43 = pkin(4) * t202 - qJ(5) * t75;
t136 = t185 * t191 + t186 * t189;
t86 = t188 * t135 + t136 * t277;
t283 = t86 / 0.2e1;
t316 = -m(6) - m(5);
t233 = qJD(1) * qJD(2);
t222 = t189 * t233;
t232 = qJDD(1) * t191;
t143 = -t222 + t232;
t315 = t143 / 0.2e1;
t312 = qJD(2) / 0.2e1;
t276 = mrSges(5,3) * t75;
t58 = -mrSges(5,2) * t183 + t276;
t61 = mrSges(6,2) * t75 + mrSges(6,3) * t183;
t259 = -t58 - t61;
t275 = mrSges(5,3) * t202;
t306 = mrSges(6,2) * t202 + t322 * t183 + t275;
t249 = qJDD(2) / 0.2e1;
t301 = t168 * pkin(4) + t167 * qJ(5);
t171 = pkin(6) * t232;
t133 = -pkin(6) * t222 + t171;
t144 = qJDD(1) * t189 + t191 * t233;
t134 = t144 * pkin(6);
t299 = t133 * t191 + t134 * t189;
t297 = 0.2e1 * t249;
t19 = t188 * t53 + t277 * t55;
t80 = qJDD(2) * pkin(2) - qJ(3) * t144 - qJD(3) * t236 - t134;
t238 = qJD(2) * t189;
t229 = pkin(6) * t238;
t237 = qJD(3) * t191;
t87 = qJ(3) * t143 + t171 + (-t229 + t237) * qJD(1);
t46 = -t185 * t87 + t186 * t80;
t92 = t143 * t185 + t144 * t186;
t21 = qJDD(2) * pkin(3) - pkin(7) * t92 + t46;
t47 = t185 * t80 + t186 * t87;
t91 = t143 * t186 - t144 * t185;
t32 = pkin(7) * t91 + t47;
t6 = -qJD(4) * t19 - t188 * t32 + t21 * t277;
t295 = -m(3) * pkin(6) + m(4) * t187 + mrSges(2,2) - mrSges(6,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t294 = m(3) * pkin(1) + m(4) * t170 + mrSges(2,1) - t321;
t12 = t183 * qJ(5) + t19;
t67 = Ifges(6,5) * t202;
t39 = t183 * Ifges(6,6) - Ifges(6,3) * t75 + t67;
t274 = Ifges(5,4) * t202;
t40 = Ifges(5,2) * t75 + t183 * Ifges(5,6) + t274;
t292 = t19 * mrSges(5,3) + t12 * mrSges(6,2) - t39 / 0.2e1 + t40 / 0.2e1 - t31 * mrSges(6,1) - t93 * mrSges(5,1);
t288 = t75 / 0.2e1;
t287 = -t75 / 0.2e1;
t285 = -t202 / 0.2e1;
t284 = t202 / 0.2e1;
t281 = -t122 / 0.2e1;
t279 = -t183 / 0.2e1;
t278 = t183 / 0.2e1;
t270 = pkin(2) * t189;
t269 = pkin(6) * t191;
t257 = mrSges(5,2) * t168;
t256 = mrSges(4,3) * t121;
t255 = mrSges(4,3) * t122;
t254 = Ifges(3,4) * t189;
t253 = Ifges(3,4) * t191;
t252 = Ifges(4,4) * t122;
t248 = qJDD(1) * pkin(1);
t218 = qJD(2) * t187;
t118 = t189 * t218 + t237;
t119 = -qJD(3) * t189 + t191 * t218;
t64 = t186 * t118 + t185 * t119;
t95 = t185 * t154 - t186 * t156;
t239 = pkin(3) * t178 + t180;
t231 = t188 * t272;
t230 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t236) * t269;
t174 = pkin(2) * t238;
t173 = pkin(2) * t236;
t226 = -t91 * mrSges(4,1) + t92 * mrSges(4,2);
t29 = qJD(4) * t75 + t188 * t91 + t277 * t92;
t30 = qJD(4) * t202 + t188 * t92 - t277 * t91;
t225 = t30 * mrSges(5,1) + t29 * mrSges(5,2);
t224 = t30 * mrSges(6,1) - t29 * mrSges(6,3);
t223 = qJD(4) * t277;
t220 = t190 * t317;
t219 = t192 * t317;
t181 = qJDD(2) + qJDD(4);
t14 = -t181 * mrSges(6,1) + t29 * mrSges(6,2);
t96 = -pkin(3) * t122 + t173;
t123 = t136 * qJD(2);
t97 = pkin(3) * t123 + t174;
t63 = -t118 * t185 + t186 * t119;
t94 = t186 * t154 + t156 * t185;
t101 = -pkin(3) * t135 - t170;
t214 = mrSges(3,1) * t189 + mrSges(3,2) * t191;
t210 = t191 * Ifges(3,2) + t254;
t209 = Ifges(3,5) * t191 - Ifges(3,6) * t189;
t105 = -qJD(4) * t231 + t169 * t223;
t124 = t135 * qJD(2);
t207 = -pkin(7) * t124 + t63;
t65 = -pkin(7) * t136 + t94;
t66 = pkin(7) * t135 + t95;
t204 = -t188 * t66 + t277 * t65;
t38 = t188 * t65 + t277 * t66;
t203 = pkin(1) * t214;
t112 = -pkin(2) * t143 + qJDD(3) - t248;
t5 = -qJD(4) * t250 + t188 * t21 + t53 * t223 + t277 * t32;
t201 = t135 * t277 - t188 * t136;
t200 = t189 * (Ifges(3,1) * t191 - t254);
t116 = t169 * t277 - t231;
t145 = -pkin(3) * t177 - t270;
t196 = m(6) * (-pkin(4) * t167 + t145) - t167 * mrSges(6,1);
t194 = t257 + (m(6) * pkin(4) - t322) * t167;
t56 = -pkin(3) * t91 + t112;
t2 = qJ(5) * t181 + qJD(5) * t183 + t5;
t3 = -t181 * pkin(4) + qJDD(5) - t6;
t193 = t6 * mrSges(5,1) - t3 * mrSges(6,1) - t5 * mrSges(5,2) + t2 * mrSges(6,3) + t308 * t30 + t309 * t29 + (Ifges(6,2) + Ifges(5,3)) * t181;
t182 = -pkin(7) + t187;
t172 = Ifges(3,4) * t235;
t153 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t235;
t140 = pkin(1) + t239;
t131 = Ifges(3,1) * t236 + Ifges(3,5) * qJD(2) + t172;
t130 = Ifges(3,6) * qJD(2) + qJD(1) * t210;
t115 = Ifges(4,4) * t121;
t114 = -pkin(4) - t116;
t113 = qJ(5) + t117;
t100 = qJD(5) + t105;
t99 = qJD(2) * mrSges(4,1) + t255;
t98 = -qJD(2) * mrSges(4,2) + t256;
t81 = -mrSges(4,1) * t121 - mrSges(4,2) * t122;
t78 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t92;
t77 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t91;
t72 = -t122 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t115;
t71 = t121 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t252;
t52 = -pkin(7) * t123 + t64;
t49 = qJD(4) * t86 + t123 * t277 + t188 * t124;
t48 = qJD(4) * t201 - t188 * t123 + t124 * t277;
t45 = -mrSges(5,1) * t75 + mrSges(5,2) * t202;
t44 = -mrSges(6,1) * t75 - mrSges(6,3) * t202;
t36 = -pkin(4) * t201 - qJ(5) * t86 + t101;
t35 = t43 + t96;
t34 = t188 * t206 + t277 * t57;
t16 = -mrSges(6,2) * t30 + mrSges(6,3) * t181;
t15 = -mrSges(5,2) * t181 - mrSges(5,3) * t30;
t13 = mrSges(5,1) * t181 - mrSges(5,3) * t29;
t10 = pkin(4) * t49 - qJ(5) * t48 - qJD(5) * t86 + t97;
t7 = pkin(4) * t30 - qJ(5) * t29 - qJD(5) * t202 + t56;
t1 = [(t201 * t5 - t6 * t86) * mrSges(5,3) + t7 * (-mrSges(6,1) * t201 - mrSges(6,3) * t86) + t56 * (-mrSges(5,1) * t201 + mrSges(5,2) * t86) + t201 * (Ifges(5,4) * t29 + Ifges(5,6) * t181) / 0.2e1 + (-t201 * t308 + t309 * t86) * t181 / 0.2e1 - t201 * (Ifges(6,5) * t29 + Ifges(6,6) * t181) / 0.2e1 + (-t201 * t310 + t311 * t86) * t29 / 0.2e1 + (t2 * t201 + t3 * t86) * mrSges(6,2) + (t316 * (t192 * t140 - t190 * t182) + (-m(6) * t301 - t294) * t192 + t295 * t190) * g(2) + ((-t182 * t316 + t295) * t192 + (-m(6) * (-t140 - t301) + m(5) * t140 + t294) * t190) * g(1) - qJDD(2) * mrSges(3,2) * t269 + Ifges(3,6) * t191 * t249 + t144 * t253 / 0.2e1 + (-t123 * t83 - t124 * t82) * mrSges(4,3) + (Ifges(4,1) * t124 - Ifges(4,4) * t123) * t281 + t121 * (Ifges(4,4) * t124 - Ifges(4,2) * t123) / 0.2e1 + t147 * (mrSges(4,1) * t123 + mrSges(4,2) * t124) - t155 * t248 + (-m(5) * t18 + m(6) * t11 + t306) * (qJD(4) * t38 + t188 * t52 - t207 * t277) + (m(5) * t19 + m(6) * t12 - t259) * (qJD(4) * t204 + t188 * t207 + t277 * t52) + (Ifges(5,4) * t288 + Ifges(6,5) * t287 + t278 * t309 + t284 * t311 - t320) * t48 - t130 * t238 / 0.2e1 - t203 * t233 + m(4) * (-t112 * t170 + t147 * t174 + t46 * t94 + t47 * t95 + t63 * t82 + t64 * t83) - t153 * t229 + (t143 * t269 + t299) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t299) + (-t112 * mrSges(4,1) + t47 * mrSges(4,3) + Ifges(4,4) * t92 + Ifges(4,2) * t91 + Ifges(4,6) * t297) * t135 + (t112 * mrSges(4,2) - t46 * mrSges(4,3) + Ifges(4,1) * t92 + Ifges(4,4) * t91 + Ifges(4,5) * t297) * t136 + t36 * t224 + t101 * t225 - t170 * t226 + t81 * t174 + (Ifges(3,1) * t144 - pkin(6) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t144) + Ifges(3,4) * t315 + t297 * Ifges(3,5)) * t189 + (-Ifges(5,2) * t288 + Ifges(6,3) * t287 + t278 * t308 + t284 * t310 - t292) * t49 + (t181 * t309 + t29 * t311) * t283 + (t209 * t312 - t230) * qJD(2) - (-m(5) * t6 + m(6) * t3 - t13 + t14) * t204 + (t200 + t191 * (-Ifges(3,2) * t189 + t253)) * t233 / 0.2e1 + m(5) * (t101 * t56 + t93 * t97) + m(6) * (t10 * t31 + t36 * t7) + t94 * t78 + t95 * t77 + t97 * t45 + t64 * t98 + t63 * t99 + (-(Ifges(6,3) + Ifges(5,2)) * t201 + 0.2e1 * t310 * t283) * t30 - t123 * t71 / 0.2e1 + t124 * t72 / 0.2e1 + (Ifges(4,5) * t124 - Ifges(4,6) * t123 + t191 * t131) * t312 + t210 * t315 - pkin(1) * (-mrSges(3,1) * t143 + mrSges(3,2) * t144) + t191 * (Ifges(3,4) * t144 + Ifges(3,2) * t143 + Ifges(3,6) * qJDD(2)) / 0.2e1 + Ifges(2,3) * qJDD(1) + t10 * t44 + (m(5) * t5 + m(6) * t2 + t15 + t16) * t38; (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t230 + (t203 - t200 / 0.2e1) * qJD(1)) * qJD(1) + (-t147 * t173 - t82 * t89 - t83 * t90 + (t185 * t47 + t186 * t46) * pkin(2)) * m(4) + t71 * t281 + t82 * t256 + t78 * t271 + t77 * t272 + t122 * (Ifges(4,1) * t121 + t252) / 0.2e1 - t83 * t255 + t259 * t34 + (Ifges(5,4) * t287 + Ifges(6,5) * t288 + t309 * t279 + t311 * t285 + t320) * t75 - t209 * t233 / 0.2e1 + (t130 / 0.2e1 + pkin(6) * t153) * t236 - t81 * t173 + (t113 * t2 + t114 * t3 - t31 * t35 + (t100 - t34) * t12 + t305 * t11) * m(6) + (t116 * t6 + t117 * t5 - t93 * t96 + (t105 - t34) * t19 - t305 * t18) * m(5) - g(1) * (t192 * t196 - t219) - g(2) * (t190 * t196 - t220) + (-Ifges(5,2) * t287 + Ifges(6,3) * t288 + t279 * t308 + t310 * t285 + t292) * t202 + (-m(4) * t180 - m(6) * (t239 + t301) - m(5) * t239 + t321) * g(3) - (-Ifges(3,2) * t236 + t131 + t172) * t235 / 0.2e1 - (Ifges(4,2) * t122 + t115 + t72) * t121 / 0.2e1 + t193 + t318 * (m(4) * t270 - m(5) * t145 + mrSges(4,1) * t177 + mrSges(5,1) * t167 + mrSges(4,2) * t178 + t214 + t257) + Ifges(4,6) * t91 + Ifges(4,5) * t92 - t96 * t45 + t305 * t306 - t90 * t98 - t89 * t99 + t100 * t61 + t105 * t58 + t113 * t16 + t114 * t14 + t116 * t13 + t117 * t15 - qJD(2) * (Ifges(4,5) * t121 + Ifges(4,6) * t122) / 0.2e1 - t133 * mrSges(3,2) - t134 * mrSges(3,1) + Ifges(3,6) * t143 + Ifges(3,5) * t144 - t147 * (-mrSges(4,1) * t122 + mrSges(4,2) * t121) - t35 * t44 + t46 * mrSges(4,1) - t47 * mrSges(4,2); -t121 * t98 - t122 * t99 - t306 * t202 + t259 * t75 + t224 + t225 + t226 + (-g(1) * t190 + g(2) * t192) * (m(4) - t316) + (-t11 * t202 - t12 * t75 + t7) * m(6) + (t18 * t202 - t19 * t75 + t56) * m(5) + (-t121 * t83 - t122 * t82 + t112) * m(4); t40 * t284 + (Ifges(6,3) * t202 + t273) * t288 + (-t306 + t275) * t19 + (t259 + t276) * t18 + t300 * g(3) + (t192 * t194 + t219) * g(1) + (t190 * t194 + t220) * g(2) + t193 - t31 * (mrSges(6,1) * t202 - mrSges(6,3) * t75) - t93 * (mrSges(5,1) * t202 + mrSges(5,2) * t75) - pkin(4) * t14 + qJ(5) * t16 - t43 * t44 + qJD(5) * t61 + (-t11 * t75 + t12 * t202) * mrSges(6,2) + (t202 * t308 + t309 * t75) * t279 + (-pkin(4) * t3 - g(3) * t301 + qJ(5) * t2 - t11 * t19 + t12 * t304 - t31 * t43) * m(6) + (-Ifges(5,2) * t202 + t307 + t68) * t287 + (t311 * t75 - t274 + t39 + t67) * t285; -t183 * t61 + t202 * t44 + (g(3) * t168 - t12 * t183 - t318 * t167 + t31 * t202 + t3) * m(6) + t14;];
tau = t1;
