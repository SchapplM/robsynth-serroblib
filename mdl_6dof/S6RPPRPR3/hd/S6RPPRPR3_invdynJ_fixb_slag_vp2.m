% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:00
% EndTime: 2019-03-09 01:44:12
% DurationCPUTime: 9.72s
% Computational Cost: add. (4713->501), mult. (8801->663), div. (0->0), fcn. (5459->14), ass. (0->223)
t302 = m(6) + m(7);
t301 = -mrSges(6,2) + mrSges(7,3);
t142 = sin(qJ(6));
t145 = cos(qJ(6));
t177 = -mrSges(7,1) * t145 + mrSges(7,2) * t142;
t300 = -m(7) * pkin(5) + t177;
t137 = qJ(1) + pkin(9);
t130 = sin(t137);
t132 = cos(t137);
t288 = -g(1) * t130 + g(2) * t132;
t136 = qJ(4) + pkin(10);
t129 = sin(t136);
t131 = cos(t136);
t143 = sin(qJ(4));
t263 = cos(qJ(4));
t165 = t143 * mrSges(5,1) + mrSges(5,2) * t263;
t298 = -mrSges(6,1) * t129 + t301 * t131 - t165;
t297 = -m(5) - m(4);
t138 = sin(pkin(10));
t230 = cos(pkin(10));
t155 = -t138 * t263 - t143 * t230;
t92 = t155 * qJD(1);
t85 = Ifges(6,4) * t92;
t181 = t230 * t263;
t215 = t143 * qJD(1);
t91 = -qJD(1) * t181 + t138 * t215;
t70 = qJD(4) * t145 + t142 * t91;
t295 = t70 * Ifges(7,6);
t290 = qJD(6) - t92;
t294 = t290 * Ifges(7,3);
t293 = t92 * Ifges(6,2);
t213 = qJD(1) * qJD(4);
t104 = qJDD(1) * t263 - t143 * t213;
t190 = qJD(4) * t263;
t105 = -qJD(1) * t190 - t143 * qJDD(1);
t64 = t104 * t230 + t138 * t105;
t34 = qJD(6) * t70 + qJDD(4) * t142 + t145 * t64;
t71 = qJD(4) * t142 - t145 * t91;
t35 = -qJD(6) * t71 + qJDD(4) * t145 - t142 * t64;
t11 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t54 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t64;
t292 = t11 - t54;
t262 = mrSges(6,3) * t91;
t291 = -qJD(4) * mrSges(6,1) - mrSges(7,1) * t70 + mrSges(7,2) * t71 - t262;
t101 = t138 * t143 - t181;
t217 = qJD(6) * t145;
t94 = t155 * qJD(4);
t163 = t101 * t217 - t142 * t94;
t218 = qJD(6) * t142;
t162 = t101 * t218 + t145 * t94;
t140 = cos(pkin(9));
t126 = -pkin(1) * t140 - pkin(2);
t118 = -pkin(7) + t126;
t98 = qJDD(1) * t118 + qJDD(3);
t99 = qJD(1) * t118 + qJD(3);
t210 = t263 * qJDD(2) + t143 * t98 + t99 * t190;
t212 = qJD(2) * qJD(4);
t42 = -t143 * t212 + t210;
t189 = t263 * qJD(2);
t74 = t143 * t99 + t189;
t43 = -qJD(4) * t74 - t143 * qJDD(2) + t263 * t98;
t289 = t42 * t143 + t43 * t263;
t231 = qJDD(4) / 0.2e1;
t287 = -t231 - qJDD(4) / 0.2e1;
t204 = mrSges(5,3) * t215;
t109 = -qJD(4) * mrSges(5,2) - t204;
t191 = qJD(1) * t263;
t183 = mrSges(5,3) * t191;
t110 = qJD(4) * mrSges(5,1) - t183;
t219 = qJD(4) * t143;
t88 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t105;
t286 = t109 * t190 - t110 * t219 + t143 * t88;
t139 = sin(pkin(9));
t121 = pkin(1) * t139 + qJ(3);
t108 = t121 * qJD(1);
t285 = -qJD(1) * t108 + t288;
t284 = mrSges(3,1) + mrSges(5,3) + mrSges(6,3) - mrSges(4,2);
t67 = t189 + (-qJ(5) * qJD(1) + t99) * t143;
t60 = t230 * t67;
t73 = -t143 * qJD(2) + t263 * t99;
t66 = -qJ(5) * t191 + t73;
t62 = qJD(4) * pkin(4) + t66;
t32 = t138 * t62 + t60;
t28 = qJD(4) * pkin(8) + t32;
t90 = pkin(4) * t215 + qJD(5) + t108;
t41 = -pkin(5) * t92 + pkin(8) * t91 + t90;
t12 = -t142 * t28 + t145 * t41;
t63 = -t138 * t104 + t105 * t230;
t214 = qJD(1) * qJD(3);
t100 = qJDD(1) * t121 + t214;
t65 = -pkin(4) * t105 + qJDD(5) + t100;
t18 = -pkin(5) * t63 - pkin(8) * t64 + t65;
t188 = t263 * qJD(5);
t26 = qJDD(4) * pkin(4) - t104 * qJ(5) - qJD(1) * t188 + t43;
t33 = qJ(5) * t105 + (-qJD(1) * qJD(5) - t212) * t143 + t210;
t10 = t138 * t26 + t230 * t33;
t8 = qJDD(4) * pkin(8) + t10;
t1 = qJD(6) * t12 + t142 * t18 + t145 * t8;
t13 = t142 * t41 + t145 * t28;
t172 = t12 * t145 + t13 * t142;
t2 = -qJD(6) * t13 - t142 * t8 + t145 * t18;
t283 = qJD(6) * t172 - t1 * t145 + t142 * t2;
t240 = t138 * t67;
t31 = t230 * t62 - t240;
t9 = -t138 * t33 + t230 * t26;
t93 = -qJD(4) * t181 + t138 * t219;
t282 = t10 * t155 + t101 * t9 - t31 * t94 + t32 * t93;
t281 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t256 = pkin(8) * t131;
t280 = mrSges(3,2) - m(7) * (pkin(5) * t129 - t256) - mrSges(4,3) + t298;
t7 = -qJDD(4) * pkin(5) - t9;
t278 = -m(6) * t9 + m(7) * t7 + t292;
t27 = -qJD(4) * pkin(5) - t31;
t277 = m(6) * t31 - m(7) * t27 - t291;
t59 = qJDD(6) - t63;
t16 = mrSges(7,1) * t59 - mrSges(7,3) * t34;
t17 = -mrSges(7,2) * t59 + mrSges(7,3) * t35;
t44 = -mrSges(7,2) * t290 + mrSges(7,3) * t70;
t45 = mrSges(7,1) * t290 - mrSges(7,3) * t71;
t276 = -m(7) * t283 - t142 * t16 + t145 * t17 - t45 * t217 - t44 * t218;
t275 = qJD(1) ^ 2;
t274 = t34 / 0.2e1;
t273 = t35 / 0.2e1;
t272 = t59 / 0.2e1;
t271 = -t70 / 0.2e1;
t270 = -t71 / 0.2e1;
t269 = t71 / 0.2e1;
t268 = -t290 / 0.2e1;
t267 = -t91 / 0.2e1;
t264 = t142 / 0.2e1;
t261 = mrSges(6,3) * t92;
t260 = Ifges(6,4) * t91;
t259 = Ifges(7,4) * t71;
t144 = sin(qJ(1));
t258 = pkin(1) * t144;
t257 = pkin(4) * t138;
t253 = t101 * t7;
t134 = t143 * pkin(4);
t146 = cos(qJ(1));
t135 = t146 * pkin(1);
t246 = Ifges(5,4) * t143;
t245 = Ifges(7,4) * t142;
t244 = Ifges(7,4) * t145;
t239 = t142 * mrSges(7,3);
t235 = t145 * mrSges(7,3);
t228 = t101 * t142;
t227 = t101 * t145;
t226 = t130 * t142;
t225 = t130 * t145;
t224 = t132 * t142;
t223 = t132 * t145;
t222 = qJ(5) - t118;
t216 = qJDD(1) * mrSges(4,2);
t115 = pkin(4) * t190 + qJD(3);
t211 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t59;
t206 = Ifges(5,4) * t263;
t23 = Ifges(7,2) * t70 + Ifges(7,6) * t290 + t259;
t201 = t23 * t264;
t69 = Ifges(7,4) * t70;
t24 = t71 * Ifges(7,1) + Ifges(7,5) * t290 + t69;
t200 = -t145 * t24 / 0.2e1;
t199 = t132 * pkin(2) + t130 * qJ(3) + t135;
t198 = t263 * t118;
t197 = t230 * pkin(4);
t192 = -t63 * mrSges(6,1) + t64 * mrSges(6,2);
t187 = t217 / 0.2e1;
t186 = t132 * qJ(3) - t258;
t185 = -t213 / 0.2e1;
t184 = pkin(4) * t191;
t107 = t121 + t134;
t180 = -g(1) * t132 - g(2) * t130;
t176 = mrSges(7,1) * t142 + mrSges(7,2) * t145;
t175 = Ifges(7,1) * t145 - t245;
t174 = -Ifges(7,2) * t142 + t244;
t173 = Ifges(7,5) * t145 - Ifges(7,6) * t142;
t171 = t12 * t142 - t13 * t145;
t48 = -pkin(5) * t155 + pkin(8) * t101 + t107;
t159 = -qJ(5) * t263 + t198;
t95 = t222 * t143;
t52 = t138 * t159 - t230 * t95;
t19 = -t142 * t52 + t145 * t48;
t20 = t142 * t48 + t145 * t52;
t170 = t108 * qJD(3) + t100 * t121;
t169 = -pkin(2) * t130 + t186;
t168 = mrSges(5,1) * t263 - mrSges(5,2) * t143;
t167 = t263 * Ifges(5,1) - t246;
t166 = -Ifges(5,2) * t143 + t206;
t164 = -Ifges(5,5) * t143 - Ifges(5,6) * t263;
t158 = t108 * t168;
t157 = t143 * (-Ifges(5,2) * t263 - t246);
t153 = (-Ifges(5,1) * t143 - t206) * t263;
t150 = t219 * t222 - t188;
t148 = (-t143 * t73 + t263 * t74) * qJD(4) + t289;
t141 = -qJ(5) - pkin(7);
t125 = -t197 - pkin(5);
t106 = qJDD(1) * t126 + qJDD(3);
t103 = t165 * qJD(1);
t97 = Ifges(5,5) * qJD(4) + qJD(1) * t167;
t96 = Ifges(5,6) * qJD(4) + qJD(1) * t166;
t87 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t104;
t79 = -qJD(4) * mrSges(6,2) + t261;
t78 = t129 * t223 - t226;
t77 = t129 * t224 + t225;
t76 = t129 * t225 + t224;
t75 = -t129 * t226 + t223;
t72 = qJD(4) * t159 - t143 * qJD(5);
t55 = -mrSges(6,1) * t92 - mrSges(6,2) * t91;
t53 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t63;
t50 = -t91 * Ifges(6,1) + Ifges(6,5) * qJD(4) + t85;
t49 = Ifges(6,6) * qJD(4) - t260 + t293;
t47 = -t91 * pkin(5) - t92 * pkin(8) + t184;
t46 = -pkin(5) * t93 - pkin(8) * t94 + t115;
t39 = t138 * t150 + t230 * t72;
t37 = t230 * t66 - t240;
t36 = t138 * t66 + t60;
t22 = t71 * Ifges(7,5) + t294 + t295;
t15 = t142 * t47 + t145 * t37;
t14 = -t142 * t37 + t145 * t47;
t6 = t34 * Ifges(7,1) + t35 * Ifges(7,4) + t59 * Ifges(7,5);
t5 = t34 * Ifges(7,4) + t35 * Ifges(7,2) + t59 * Ifges(7,6);
t4 = -qJD(6) * t20 - t142 * t39 + t145 * t46;
t3 = qJD(6) * t19 + t142 * t46 + t145 * t39;
t21 = [(m(3) * t258 - m(4) * t169 - m(5) * t186 + mrSges(2,1) * t144 - t78 * mrSges(7,1) + mrSges(2,2) * t146 + t77 * mrSges(7,2) - t302 * (t130 * t141 + t132 * t134 + t169) + (-m(5) * (-pkin(2) - pkin(7)) + t284) * t130 + t280 * t132) * g(1) + (-m(3) * t135 - mrSges(2,1) * t146 - t76 * mrSges(7,1) + mrSges(2,2) * t144 - t75 * mrSges(7,2) + t297 * t199 - t302 * (t130 * t134 - t132 * t141 + t199) + (-m(5) * pkin(7) - t284) * t132 + t280 * t130) * g(2) + m(7) * (t1 * t20 + t12 * t4 + t13 * t3 + t19 * t2) + m(6) * (t10 * t52 + t107 * t65 + t115 * t90 + t32 * t39) + (t121 * mrSges(4,3) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (0.2e1 * mrSges(3,1) * t140 - 0.2e1 * mrSges(3,2) * t139 + m(3) * (t139 ^ 2 + t140 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t290 * (Ifges(7,5) * t162 + Ifges(7,6) * t163) / 0.2e1 - (t211 / 0.2e1 - Ifges(6,4) * t64 - Ifges(6,2) * t63 + t65 * mrSges(6,1) + Ifges(7,3) * t272 + Ifges(7,6) * t273 + Ifges(7,5) * t274 + t287 * Ifges(6,6) + t281) * t155 + t100 * t165 + (Ifges(5,5) * t263 - Ifges(5,6) * t143) * t231 + t263 * (Ifges(5,1) * t104 + Ifges(5,4) * t105 + Ifges(5,5) * qJDD(4)) / 0.2e1 - t176 * t253 + m(4) * (t106 * t126 + t170) + t282 * mrSges(6,3) - (t200 + t201 - t50 / 0.2e1 - t90 * mrSges(6,2) - t85 / 0.2e1 - Ifges(6,1) * t267) * t94 + (-t12 * mrSges(7,1) + t13 * mrSges(7,2) - t294 / 0.2e1 - t295 / 0.2e1 + t49 / 0.2e1 - t22 / 0.2e1 - t90 * mrSges(6,1) + t293 / 0.2e1 + Ifges(6,4) * t267 - Ifges(7,5) * t269) * t93 + (t158 + Ifges(6,5) * t94 / 0.2e1 + Ifges(6,6) * t93 / 0.2e1 + t164 * qJD(4) / 0.2e1) * qJD(4) + (-t190 * t74 + t219 * t73 - t289) * mrSges(5,3) + (-t65 * mrSges(6,2) - Ifges(6,1) * t64 - Ifges(6,4) * t63 + Ifges(6,5) * t287 - t173 * t272 - t174 * t273 - t175 * t274 + t23 * t187) * t101 + (m(5) * t148 + t286) * t118 - t6 * t227 / 0.2e1 - t97 * t219 / 0.2e1 + t153 * t213 / 0.2e1 + (qJD(6) * t24 + t5) * t228 / 0.2e1 - t96 * t190 / 0.2e1 + t107 * t192 + t70 * (Ifges(7,4) * t162 + Ifges(7,2) * t163) / 0.2e1 + t105 * t166 / 0.2e1 + t104 * t167 / 0.2e1 + t27 * (-mrSges(7,1) * t163 + mrSges(7,2) * t162) + (Ifges(7,1) * t162 + Ifges(7,4) * t163) * t269 + m(5) * t170 + t3 * t44 + t4 * t45 + t19 * t16 + t20 * t17 + (t1 * t228 - t12 * t162 + t13 * t163 + t2 * t227) * mrSges(7,3) + t157 * t185 - t277 * (t138 * t72 - t150 * t230) + t278 * (-t138 * t95 - t159 * t230) + t126 * t216 + (t214 + t100) * mrSges(4,3) + t87 * t198 + t52 * t53 + t39 * t79 + qJD(3) * t103 + t106 * mrSges(4,2) + t115 * t55 + t121 * (-mrSges(5,1) * t105 + mrSges(5,2) * t104) - t143 * (Ifges(5,4) * t104 + Ifges(5,2) * t105 + Ifges(5,6) * qJDD(4)) / 0.2e1; m(5) * (t42 * t263 - t43 * t143 + (-t143 * t74 - t263 * t73) * qJD(4)) + t263 * t88 + m(7) * (t101 * t283 - t171 * t94) - t17 * t227 - t109 * t219 - t110 * t190 + m(6) * (-t10 * t101 + t32 * t94) + t16 * t228 + t94 * t79 - t101 * t53 - t143 * t87 + t277 * t93 + t163 * t45 + t162 * t44 - t278 * t155 + (m(4) + m(3)) * qJDD(2) + (-m(3) - t302 + t297) * g(3); -t275 * mrSges(4,3) + t263 * t87 + t216 - t291 * t94 + t292 * t101 + (t142 * t45 - t145 * t44 - t79) * t93 + (t171 * t93 - t27 * t94 + t253 + t288) * m(7) + (-t282 + t288) * m(6) + (t148 + t285) * m(5) + (t106 + t285) * m(4) - (t276 + t53) * t155 + (-m(6) * t90 - m(7) * t172 - t142 * t44 - t145 * t45 - t103 - t55) * qJD(1) + t286; (t173 * t290 + t174 * t70 + t175 * t71) * qJD(6) / 0.2e1 - t32 * t262 + (-m(7) * (-t134 + t256) - t300 * t129 + m(6) * t134 - t298) * g(3) - t55 * t184 + t7 * t177 + (-t12 * t14 + t125 * t7 - t13 * t15 - t27 * t36) * m(7) - t291 * t36 + (t168 + t302 * t263 * pkin(4) + (mrSges(6,1) - t300) * t131 + (m(7) * pkin(8) + t301) * t129) * t288 - t12 * (-mrSges(7,1) * t91 - t235 * t92) - t23 * t218 / 0.2e1 + t97 * t215 / 0.2e1 + t290 * t27 * t176 + (t183 + t110) * t74 - t13 * (mrSges(7,2) * t91 - t239 * t92) - t2 * t239 - (Ifges(6,2) * t91 + t50 + t85) * t92 / 0.2e1 + (Ifges(6,1) * t92 + t22 + t260) * t91 / 0.2e1 + t96 * t191 / 0.2e1 + (-t204 - t109) * t73 - qJD(1) * t158 - t42 * mrSges(5,2) + t43 * mrSges(5,1) - t15 * t44 - t14 * t45 + t9 * mrSges(6,1) - t10 * mrSges(6,2) + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) + (t157 / 0.2e1 - t153 / 0.2e1) * t275 + (-t12 * t217 - t13 * t218) * mrSges(7,3) + t92 * t200 + t92 * t201 + t24 * t187 + t164 * t185 + t276 * (pkin(8) + t257) + (-t184 * t90 + t31 * t36 - t32 * t37 + (t10 * t138 + t230 * t9) * pkin(4)) * m(6) + t1 * t235 + t54 * t197 + Ifges(6,6) * t63 + Ifges(6,5) * t64 - t37 * t79 - qJD(4) * (Ifges(6,5) * t92 + Ifges(6,6) * t91) / 0.2e1 - t90 * (-mrSges(6,1) * t91 + mrSges(6,2) * t92) + Ifges(5,5) * t104 + Ifges(5,6) * t105 + t125 * t11 + t145 * t5 / 0.2e1 + t53 * t257 + t31 * t261 + t6 * t264 + t49 * t267 + (-Ifges(7,3) * t91 + t173 * t92) * t268 + (-Ifges(7,5) * t91 + t175 * t92) * t270 + (-Ifges(7,6) * t91 + t174 * t92) * t271 + (Ifges(7,5) * t142 + Ifges(7,6) * t145) * t272 + (Ifges(7,2) * t145 + t245) * t273 + (Ifges(7,1) * t142 + t244) * t274; -t92 * t79 + t291 * t91 + (t290 * t44 + t16) * t145 + (-t290 * t45 + t17) * t142 + t192 + (t1 * t142 + t2 * t145 - t171 * t290 + t27 * t91 + t180) * m(7) + (-t31 * t91 - t32 * t92 + t180 + t65) * m(6); -t27 * (mrSges(7,1) * t71 + mrSges(7,2) * t70) + (Ifges(7,1) * t70 - t259) * t270 + t23 * t269 + (Ifges(7,5) * t70 - Ifges(7,6) * t71) * t268 - t12 * t44 + t13 * t45 - g(1) * (mrSges(7,1) * t75 - mrSges(7,2) * t76) - g(2) * (mrSges(7,1) * t77 + mrSges(7,2) * t78) + g(3) * t176 * t131 + (t12 * t70 + t13 * t71) * mrSges(7,3) + t211 + (-Ifges(7,2) * t71 + t24 + t69) * t271 + t281;];
tau  = t21;
