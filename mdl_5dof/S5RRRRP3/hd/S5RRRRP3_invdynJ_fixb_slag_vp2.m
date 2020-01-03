% Calculate vector of inverse dynamics joint torques for
% S5RRRRP3
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:07
% EndTime: 2019-12-31 21:49:14
% DurationCPUTime: 3.32s
% Computational Cost: add. (3122->300), mult. (4468->377), div. (0->0), fcn. (2165->12), ass. (0->150)
t289 = mrSges(5,1) + mrSges(6,1);
t284 = Ifges(6,4) + Ifges(5,5);
t283 = Ifges(6,6) - Ifges(5,6);
t132 = qJD(1) + qJD(2);
t126 = qJD(3) + t132;
t136 = sin(qJ(4));
t140 = cos(qJ(4));
t213 = t126 * t140;
t195 = mrSges(6,2) * t213;
t86 = qJD(4) * mrSges(6,3) + t195;
t271 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t213 + t86;
t214 = t126 * t136;
t196 = mrSges(6,2) * t214;
t272 = -mrSges(5,3) * t214 + t289 * qJD(4) - t196;
t137 = sin(qJ(3));
t141 = cos(qJ(3));
t138 = sin(qJ(2));
t233 = pkin(1) * qJD(1);
t200 = t138 * t233;
t142 = cos(qJ(2));
t209 = qJD(1) * t142;
t199 = pkin(1) * t209;
t89 = pkin(2) * t132 + t199;
t48 = t137 * t89 + t141 * t200;
t37 = pkin(8) * t126 + t48;
t228 = t136 * t37;
t28 = -qJD(4) * pkin(4) + qJD(5) + t228;
t223 = t140 * t37;
t29 = qJD(4) * qJ(5) + t223;
t288 = -t136 * t272 + t140 * t271 + m(6) * (t136 * t28 + t140 * t29) - t126 * mrSges(4,2);
t158 = -t136 * t29 + t140 * t28;
t131 = qJDD(1) + qJDD(2);
t125 = qJDD(3) + t131;
t207 = qJD(3) * t141;
t250 = pkin(1) * t142;
t81 = -qJD(2) * t200 + qJDD(1) * t250;
t281 = pkin(2) * t131 - qJD(3) * t200 + t81;
t82 = (qJD(2) * t209 + qJDD(1) * t138) * pkin(1);
t15 = t281 * t137 + t141 * t82 + t89 * t207;
t12 = pkin(8) * t125 + t15;
t11 = t140 * t12;
t205 = qJD(4) * t136;
t8 = -t205 * t37 + t11;
t204 = qJD(4) * t140;
t9 = -t12 * t136 - t204 * t37;
t171 = -t136 * t9 + t140 * t8;
t5 = qJDD(4) * qJ(5) + t11 + (qJD(5) - t228) * qJD(4);
t6 = -qJDD(4) * pkin(4) + qJDD(5) - t9;
t278 = t136 * t6 + t140 * t5;
t64 = t125 * t136 + t126 * t204;
t45 = -qJDD(4) * mrSges(6,1) + t64 * mrSges(6,2);
t63 = -t140 * t125 + t126 * t205;
t46 = -mrSges(6,2) * t63 + qJDD(4) * mrSges(6,3);
t276 = -t204 * t272 - t205 * t271 + m(5) * t171 + m(6) * (qJD(4) * t158 + t278) + (-qJDD(4) * mrSges(5,2) - mrSges(5,3) * t63 + t46) * t140 + t136 * (-qJDD(4) * mrSges(5,1) + mrSges(5,3) * t64 + t45);
t166 = t140 * mrSges(6,1) + t136 * mrSges(6,3);
t98 = -t140 * mrSges(5,1) + mrSges(5,2) * t136;
t287 = t98 - t166;
t285 = -mrSges(6,2) - mrSges(5,3) + mrSges(4,2);
t208 = qJD(3) * t137;
t16 = -t137 * t82 + t281 * t141 - t89 * t208;
t13 = -pkin(3) * t125 - t16;
t265 = m(5) * t13 + mrSges(5,1) * t63 + mrSges(5,2) * t64;
t100 = Ifges(5,4) * t213;
t234 = Ifges(6,5) * t140;
t164 = t136 * Ifges(6,1) - t234;
t282 = Ifges(5,1) * t214 + t284 * qJD(4) + t126 * t164 + t100;
t165 = t136 * mrSges(6,1) - t140 * mrSges(6,3);
t167 = mrSges(5,1) * t136 + mrSges(5,2) * t140;
t103 = t137 * t200;
t47 = t141 * t89 - t103;
t160 = pkin(4) * t140 + qJ(5) * t136;
t90 = -pkin(3) - t160;
t23 = t126 * t90 - t47;
t36 = -pkin(3) * t126 - t47;
t280 = t23 * t165 + t36 * t167;
t279 = t283 * t136 + t284 * t140;
t273 = m(4) * t47;
t270 = (-mrSges(4,1) + t98) * t126;
t252 = pkin(1) * t138;
t269 = (mrSges(3,1) * t252 + mrSges(3,2) * t250) * t132;
t135 = qJ(1) + qJ(2);
t129 = qJ(3) + t135;
t117 = sin(t129);
t118 = cos(t129);
t268 = g(1) * t118 + g(2) * t117;
t267 = t285 * t118 + (-m(6) * t90 + mrSges(4,1) - t287) * t117;
t217 = t118 * t140;
t218 = t118 * t136;
t266 = -t118 * mrSges(4,1) + (mrSges(5,2) - mrSges(6,3)) * t218 - t289 * t217 + t285 * t117;
t127 = sin(t135);
t128 = cos(t135);
t262 = mrSges(3,1) * t127 + mrSges(3,2) * t128 + t267;
t261 = -t128 * mrSges(3,1) + t127 * mrSges(3,2) + t266;
t260 = m(5) * t36 + t270;
t60 = t166 * t126;
t259 = -m(6) * t23 - t260 + t60;
t181 = (t136 ^ 2 + t140 ^ 2) * t37;
t257 = -m(5) * t181 - t288;
t256 = m(4) * t48 - t257;
t253 = t136 / 0.2e1;
t139 = sin(qJ(1));
t251 = pkin(1) * t139;
t249 = pkin(2) * t127;
t116 = pkin(2) * t128;
t248 = pkin(2) * t137;
t247 = pkin(2) * t141;
t143 = cos(qJ(1));
t130 = t143 * pkin(1);
t237 = Ifges(5,4) * t136;
t236 = Ifges(5,4) * t140;
t235 = Ifges(6,5) * t136;
t232 = t125 * mrSges(4,1);
t231 = t125 * mrSges(4,2);
t212 = t137 * t138;
t211 = t138 * t141;
t210 = t118 * pkin(3) + t117 * pkin(8);
t121 = pkin(2) + t250;
t74 = pkin(1) * t211 + t137 * t121;
t206 = qJD(4) * t126;
t203 = qJD(5) * t136;
t198 = pkin(2) * t208;
t190 = t116 + t210;
t185 = -t206 / 0.2e1;
t110 = t118 * pkin(8);
t182 = -pkin(3) * t117 + t110;
t73 = -pkin(1) * t212 + t121 * t141;
t176 = pkin(4) * t217 + qJ(5) * t218 + t210;
t172 = -t249 - t251;
t170 = t116 + t176;
t163 = Ifges(5,2) * t140 + t237;
t159 = pkin(4) * t136 - qJ(5) * t140;
t156 = t137 * t142 + t211;
t153 = t182 - t249;
t150 = t136 * (Ifges(5,1) * t140 - t237);
t149 = t140 * (Ifges(6,3) * t136 + t234);
t70 = pkin(4) * t205 - qJ(5) * t204 - t203;
t31 = t121 * t208 + (qJD(2) * t156 + t138 * t207) * pkin(1);
t2 = pkin(4) * t63 - qJ(5) * t64 - t126 * t203 + t13;
t99 = Ifges(6,5) * t214;
t50 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t213 + t99;
t51 = Ifges(5,6) * qJD(4) + t126 * t163;
t145 = t140 * (Ifges(5,4) * t64 + Ifges(5,6) * qJDD(4)) / 0.2e1 - t140 * (Ifges(6,5) * t64 + Ifges(6,6) * qJDD(4)) / 0.2e1 + Ifges(4,3) * t125 + t13 * t98 - t15 * mrSges(4,2) + t16 * mrSges(4,1) - t51 * t205 / 0.2e1 + t150 * t206 / 0.2e1 - t2 * t166 + t149 * t185 + (t136 * Ifges(5,1) + t164 + t236) * t64 / 0.2e1 + ((Ifges(5,1) + Ifges(6,1)) * t64 + t284 * qJDD(4)) * t253 + (t284 * t136 - t283 * t140) * qJDD(4) / 0.2e1 + (t126 * (Ifges(6,1) * t140 + t235) + t50) * t205 / 0.2e1 + t171 * mrSges(5,3) + (t235 / 0.2e1 - t163 / 0.2e1 + (Ifges(6,5) - Ifges(5,4)) * t253 + (-Ifges(5,2) / 0.2e1 - Ifges(6,3)) * t140) * t63 + (t126 * (-Ifges(5,2) * t136 + t236) + t282) * t204 / 0.2e1 + (t204 * t28 - t205 * t29 + t278) * mrSges(6,2) + (t280 + t279 * qJD(4) / 0.2e1) * qJD(4);
t144 = t81 * mrSges(3,1) - t82 * mrSges(3,2) + Ifges(3,3) * t131 + t145;
t80 = t90 - t247;
t62 = t159 * t126;
t49 = t70 + t198;
t42 = -t73 + t90;
t26 = mrSges(6,1) * t63 - mrSges(6,3) * t64;
t18 = t31 + t70;
t1 = [m(4) * (t15 * t74 + t16 * t73) + t144 + m(6) * (t18 * t23 + t2 * t42) + t73 * t232 - t74 * t231 - t18 * t60 + t42 * t26 + Ifges(2,3) * qJDD(1) + m(3) * (t138 * t82 + t142 * t81) * pkin(1) + t265 * (-pkin(3) - t73) + (t260 - t273) * t31 + (mrSges(3,1) * t250 - mrSges(3,2) * t252) * t131 - t269 * qJD(2) + t256 * (t121 * t207 + (-t138 * t208 + (t141 * t142 - t212) * qJD(2)) * pkin(1)) + (-m(6) * (t130 + t170) - m(5) * (t130 + t190) - m(4) * (t116 + t130) - mrSges(2,1) * t143 + t139 * mrSges(2,2) - m(3) * t130 + t261) * g(2) + (t139 * mrSges(2,1) + mrSges(2,2) * t143 + m(3) * t251 - m(5) * (t153 - t251) - m(6) * (t110 + t172) - m(4) * t172 + t262) * g(1) + t276 * (pkin(8) + t74); t144 + m(6) * (t2 * t80 + t23 * t49) + t80 * t26 - t231 * t248 + t232 * t247 - t49 * t60 + t270 * t198 + t269 * qJD(1) + (-m(4) * t116 - m(5) * t190 - m(6) * t170 + t261) * g(2) + (m(4) * t249 - m(6) * (t110 - t249) - m(5) * t153 + t262) * g(1) + (t259 + t273) * t156 * t233 - t256 * (t141 * t199 - t103) + t276 * (pkin(8) + t248) + t265 * (-pkin(3) - t247) + (m(5) * (t137 * t36 + t141 * t181) * qJD(3) + m(4) * (t137 * t15 + t141 * t16 + (-t137 * t47 + t141 * t48) * qJD(3)) + t288 * t207) * pkin(2); t145 + t90 * t26 - t70 * t60 + m(6) * (t2 * t90 + t23 * t70) - t265 * pkin(3) + (-m(5) * t210 - m(6) * t176 + t266) * g(2) + (-m(5) * t182 - m(6) * t110 + t267) * g(1) + t259 * t48 + t257 * t47 + t276 * pkin(8); t271 * t228 + t272 * t223 + t283 * t63 + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + t279 * t185 - t28 * t195 - (-Ifges(5,2) * t214 + t100 + t282) * t213 / 0.2e1 - (Ifges(6,1) * t213 + t50 + t99) * t214 / 0.2e1 + qJD(5) * t86 + t284 * t64 + t29 * t196 + t51 * t214 / 0.2e1 + t287 * g(3) + t62 * t60 - pkin(4) * t45 + qJ(5) * t46 + t5 * mrSges(6,3) - t6 * mrSges(6,1) - t8 * mrSges(5,2) + t9 * mrSges(5,1) + (-pkin(4) * t6 - g(3) * t160 + qJ(5) * t5 + t29 * qJD(5) - t158 * t37 - t23 * t62) * m(6) + (m(6) * t159 + t165 + t167) * t268 + ((-t150 / 0.2e1 + t149 / 0.2e1) * t126 - t280) * t126; -t60 * t214 - qJD(4) * t86 + (g(3) * t140 - t29 * qJD(4) - t136 * t268 + t23 * t214 + t6) * m(6) + t45;];
tau = t1;
