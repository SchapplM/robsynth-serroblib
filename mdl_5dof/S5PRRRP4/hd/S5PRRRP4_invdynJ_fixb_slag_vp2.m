% Calculate vector of inverse dynamics joint torques for
% S5PRRRP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:37
% EndTime: 2019-12-05 16:45:45
% DurationCPUTime: 3.57s
% Computational Cost: add. (2009->295), mult. (3422->384), div. (0->0), fcn. (2038->10), ass. (0->148)
t272 = mrSges(5,1) + mrSges(6,1);
t128 = qJ(2) + qJ(3);
t122 = cos(t128);
t121 = sin(t128);
t131 = sin(qJ(4));
t214 = t121 * t131;
t290 = -mrSges(6,2) - mrSges(5,3);
t291 = -m(5) - m(6);
t293 = -mrSges(5,2) * t214 + t122 * (pkin(7) * t291 + t290);
t285 = Ifges(6,4) + Ifges(5,5);
t284 = Ifges(6,6) - Ifges(5,6);
t125 = qJD(2) + qJD(3);
t134 = cos(qJ(4));
t208 = t125 * t134;
t189 = mrSges(6,2) * t208;
t89 = qJD(4) * mrSges(6,3) + t189;
t235 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t208 + t89;
t136 = cos(qJ(2));
t201 = qJD(1) * t136;
t103 = qJD(2) * pkin(2) + t201;
t132 = sin(qJ(3));
t135 = cos(qJ(3));
t133 = sin(qJ(2));
t202 = qJD(1) * t133;
t51 = t103 * t132 + t135 * t202;
t40 = pkin(7) * t125 + t51;
t224 = t131 * t40;
t27 = -qJD(4) * pkin(4) + qJD(5) + t224;
t209 = t125 * t131;
t190 = mrSges(6,2) * t209;
t270 = -mrSges(5,3) * t209 + t272 * qJD(4) - t190;
t220 = t134 * t40;
t28 = qJD(4) * qJ(5) + t220;
t286 = -t125 * mrSges(4,2) - t131 * t270 + t235 * t134 + m(6) * (t131 * t27 + t134 * t28);
t271 = mrSges(5,2) - mrSges(6,3);
t124 = qJDD(2) + qJDD(3);
t199 = qJD(3) * t135;
t194 = qJD(1) * qJD(2);
t174 = t133 * t194;
t87 = t136 * qJDD(1) - t174;
t282 = qJDD(2) * pkin(2) - qJD(3) * t202 + t87;
t173 = t136 * t194;
t88 = qJDD(1) * t133 + t173;
t15 = t103 * t199 + t282 * t132 + t135 * t88;
t12 = pkin(7) * t124 + t15;
t11 = t134 * t12;
t5 = qJDD(4) * qJ(5) + t11 + (qJD(5) - t224) * qJD(4);
t196 = qJD(4) * t134;
t9 = -t12 * t131 - t196 * t40;
t6 = -qJDD(4) * pkin(4) + qJDD(5) - t9;
t164 = t131 * t6 + t134 * t5;
t197 = qJD(4) * t131;
t289 = t27 * t196 - t28 * t197 + t164;
t200 = qJD(3) * t132;
t16 = -t103 * t200 - t132 * t88 + t282 * t135;
t13 = -pkin(3) * t124 - t16;
t66 = -t134 * t124 + t125 * t197;
t67 = t124 * t131 + t125 * t196;
t26 = mrSges(5,1) * t66 + mrSges(5,2) * t67;
t287 = m(5) * t13 + t26;
t105 = Ifges(5,4) * t208;
t230 = Ifges(6,5) * t134;
t156 = t131 * Ifges(6,1) - t230;
t283 = Ifges(5,1) * t209 + t285 * qJD(4) + t125 * t156 + t105;
t157 = t131 * mrSges(6,1) - t134 * mrSges(6,3);
t159 = mrSges(5,1) * t131 + mrSges(5,2) * t134;
t106 = t132 * t202;
t50 = t103 * t135 - t106;
t90 = -pkin(4) * t134 - qJ(5) * t131 - pkin(3);
t22 = t125 * t90 - t50;
t39 = -pkin(3) * t125 - t50;
t281 = t22 * t157 + t39 * t159;
t280 = t284 * t131 + t285 * t134;
t158 = -t134 * mrSges(6,1) - t131 * mrSges(6,3);
t234 = mrSges(5,1) * t134;
t278 = (-m(6) * t90 - t158 + t234) * t121;
t8 = -t197 * t40 + t11;
t163 = -t131 * t9 + t134 * t8;
t44 = -qJDD(4) * mrSges(6,1) + t67 * mrSges(6,2);
t45 = -mrSges(6,2) * t66 + qJDD(4) * mrSges(6,3);
t276 = m(5) * t163 + (-qJDD(4) * mrSges(5,1) + mrSges(5,3) * t67 + t44) * t131 + (-qJDD(4) * mrSges(5,2) - mrSges(5,3) * t66 + t45) * t134;
t151 = -t131 * t28 + t134 * t27;
t275 = m(6) * (qJD(4) * t151 + t164) - t197 * t235 - t196 * t270 + t276;
t274 = m(4) * t51;
t160 = mrSges(5,2) * t131 - t234;
t269 = (-mrSges(4,1) + t160) * t125;
t249 = pkin(3) * t121;
t251 = pkin(2) * t133;
t268 = m(6) * t251 - m(5) * (-t249 - t251) + t278;
t129 = sin(pkin(8));
t130 = cos(pkin(8));
t267 = g(1) * t130 + g(2) * t129;
t63 = t158 * t125;
t266 = t63 + t269;
t210 = t122 * t134;
t211 = t122 * t131;
t265 = -t122 * mrSges(4,1) + t271 * t211 - t272 * t210 + (mrSges(4,2) + t290) * t121;
t169 = (t131 ^ 2 + t134 ^ 2) * t40;
t260 = t293 * t130;
t259 = t293 * t129;
t258 = m(5) * t249 + t278;
t257 = -m(5) * t39 - m(6) * t22 - t266;
t256 = -m(5) * t169 - t286;
t253 = t131 / 0.2e1;
t252 = pkin(2) * t132;
t250 = pkin(2) * t135;
t244 = g(3) * t121;
t123 = t136 * pkin(2);
t233 = Ifges(5,4) * t131;
t232 = Ifges(5,4) * t134;
t231 = Ifges(6,5) * t131;
t229 = t124 * mrSges(4,1);
t228 = t124 * mrSges(4,2);
t207 = t129 * t131;
t206 = t129 * t134;
t205 = t130 * t131;
t204 = t130 * t134;
t203 = t122 * pkin(3) + t121 * pkin(7);
t198 = qJD(4) * t125;
t195 = qJD(5) * t131;
t188 = pkin(2) * t200;
t172 = -t198 / 0.2e1;
t165 = pkin(4) * t210 + qJ(5) * t211 + t203;
t161 = mrSges(4,1) * t121 + mrSges(4,2) * t122;
t155 = t134 * Ifges(5,2) + t233;
t152 = pkin(4) * t131 - qJ(5) * t134;
t149 = t132 * t136 + t133 * t135;
t80 = t132 * t133 - t135 * t136;
t145 = t131 * (Ifges(5,1) * t134 - t233);
t144 = t134 * (Ifges(6,3) * t131 + t230);
t62 = pkin(4) * t197 - qJ(5) * t196 - t195;
t3 = pkin(4) * t66 - qJ(5) * t67 - t125 * t195 + t13;
t104 = Ifges(6,5) * t209;
t52 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t208 + t104;
t53 = Ifges(5,6) * qJD(4) + t125 * t155;
t138 = t134 * (Ifges(5,4) * t67 + Ifges(5,6) * qJDD(4)) / 0.2e1 - t134 * (Ifges(6,5) * t67 + Ifges(6,6) * qJDD(4)) / 0.2e1 + Ifges(4,3) * t124 - t15 * mrSges(4,2) + t16 * mrSges(4,1) + t144 * t172 + t3 * t158 + t13 * t160 - t53 * t197 / 0.2e1 + t145 * t198 / 0.2e1 + (t131 * Ifges(5,1) + t156 + t232) * t67 / 0.2e1 + ((Ifges(5,1) + Ifges(6,1)) * t67 + t285 * qJDD(4)) * t253 + (t285 * t131 - t284 * t134) * qJDD(4) / 0.2e1 + (t125 * (Ifges(6,1) * t134 + t231) + t52) * t197 / 0.2e1 + t163 * mrSges(5,3) + (-t155 / 0.2e1 + t231 / 0.2e1 + (Ifges(6,5) - Ifges(5,4)) * t253 + (-Ifges(5,2) / 0.2e1 - Ifges(6,3)) * t134) * t66 + (t125 * (-Ifges(5,2) * t131 + t232) + t283) * t196 / 0.2e1 + t289 * mrSges(6,2) + (t281 + t280 * qJD(4) / 0.2e1) * qJD(4);
t137 = qJD(2) ^ 2;
t75 = t90 - t250;
t65 = t152 * t125;
t60 = t122 * t204 + t207;
t59 = t122 * t205 - t206;
t58 = t122 * t206 - t205;
t57 = t122 * t207 + t204;
t41 = t62 + t188;
t30 = t125 * t149;
t29 = t125 * t80;
t25 = mrSges(6,1) * t66 - mrSges(6,3) * t67;
t1 = [m(2) * qJDD(1) + (-qJDD(2) * t133 - t136 * t137) * mrSges(3,2) + (qJDD(2) * t136 - t133 * t137) * mrSges(3,1) + (t25 + t26 - t229) * t80 + t266 * t30 + (-m(2) - m(3) - m(4) + t291) * g(3) + m(5) * (t13 * t80 - t169 * t29 + t30 * t39) + m(6) * (t22 * t30 + t3 * t80) + m(4) * (-t16 * t80 - t30 * t50) + m(3) * (t133 * t88 + t136 * t87) + (-t228 + (-t131 * t235 - t134 * t270) * qJD(4) + m(6) * t289 + m(4) * t15 + t276) * t149 - (t274 + t286) * t29; (t173 - t88) * mrSges(3,2) + t229 * t250 + (t174 + t87) * mrSges(3,1) + (m(4) * t50 + t257) * t149 * qJD(1) + t138 + t75 * t25 + t41 * t63 + (t256 - t274) * (t135 * t201 - t106) + t269 * t188 + (t129 * t268 + t259) * g(2) + (t130 * t268 + t260) * g(1) + (-mrSges(3,1) * t136 + mrSges(3,2) * t133 - m(6) * (t123 + t165) - m(5) * (t123 + t203) - m(4) * t123 + t265) * g(3) + m(6) * (t22 * t41 + t3 * t75) + Ifges(3,3) * qJDD(2) - t228 * t252 + (m(4) * t251 + mrSges(3,1) * t133 + mrSges(3,2) * t136 + t161) * t267 + t275 * (pkin(7) + t252) + t287 * (-pkin(3) - t250) + (m(4) * (t132 * t15 + t135 * t16 + (-t132 * t50 + t135 * t51) * qJD(3)) + m(5) * (t132 * t39 + t135 * t169) * qJD(3) + t286 * t199) * pkin(2); m(6) * (t22 * t62 + t3 * t90) + t138 + t90 * t25 + t62 * t63 + t267 * t161 + (t129 * t258 + t259) * g(2) + (t130 * t258 + t260) * g(1) + (-m(5) * t203 - m(6) * t165 + t265) * g(3) + t257 * t51 + t256 * t50 - t287 * pkin(3) + t275 * pkin(7); (t159 + t157) * t244 + t28 * t190 + t284 * t66 + qJD(5) * t89 - t65 * t63 - pkin(4) * t44 + qJ(5) * t45 + (t271 * t60 + t272 * t59) * g(1) + (t271 * t58 + t272 * t57) * g(2) + t235 * t224 + t270 * t220 + t5 * mrSges(6,3) - t6 * mrSges(6,1) - t8 * mrSges(5,2) + t9 * mrSges(5,1) - (-Ifges(5,2) * t209 + t105 + t283) * t208 / 0.2e1 - (Ifges(6,1) * t208 + t104 + t52) * t209 / 0.2e1 + (-g(2) * (-pkin(4) * t57 + qJ(5) * t58) - g(1) * (-pkin(4) * t59 + qJ(5) * t60) - t151 * t40 - t22 * t65 - pkin(4) * t6 + qJ(5) * t5 + qJD(5) * t28 + t152 * t244) * m(6) + t285 * t67 - t27 * t189 + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + t53 * t209 / 0.2e1 + t280 * t172 + ((t144 / 0.2e1 - t145 / 0.2e1) * t125 - t281) * t125; t63 * t209 - qJD(4) * t89 + (-g(1) * t59 - g(2) * t57 - g(3) * t214 - t28 * qJD(4) + t209 * t22 + t6) * m(6) + t44;];
tau = t1;
