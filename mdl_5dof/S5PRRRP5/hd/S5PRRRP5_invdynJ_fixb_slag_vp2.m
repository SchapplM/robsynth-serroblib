% Calculate vector of inverse dynamics joint torques for
% S5PRRRP5
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:46
% EndTime: 2019-12-05 16:48:05
% DurationCPUTime: 7.89s
% Computational Cost: add. (2538->378), mult. (5674->496), div. (0->0), fcn. (3712->10), ass. (0->176)
t293 = mrSges(5,1) + mrSges(6,1);
t291 = mrSges(5,2) + mrSges(6,2);
t279 = Ifges(5,4) + Ifges(6,4);
t280 = Ifges(5,1) + Ifges(6,1);
t278 = Ifges(5,5) + Ifges(6,5);
t277 = Ifges(5,2) + Ifges(6,2);
t276 = Ifges(5,6) + Ifges(6,6);
t160 = sin(qJ(4));
t161 = sin(qJ(3));
t163 = cos(qJ(4));
t164 = cos(qJ(3));
t183 = t160 * t161 - t163 * t164;
t111 = t183 * qJD(2);
t290 = t279 * t111;
t121 = t160 * t164 + t161 * t163;
t112 = t121 * qJD(2);
t289 = t279 * t112;
t157 = qJ(3) + qJ(4);
t151 = cos(t157);
t243 = pkin(3) * t164;
t134 = pkin(4) * t151 + t243;
t145 = pkin(2) + t243;
t150 = sin(t157);
t189 = -mrSges(4,1) * t164 + mrSges(4,2) * t161;
t288 = mrSges(3,1) + m(6) * (pkin(2) + t134) + m(5) * t145 + m(4) * pkin(2) - t189 + t293 * t151 - t291 * t150;
t166 = -pkin(7) - pkin(6);
t287 = mrSges(3,2) + m(6) * (-qJ(5) + t166) - mrSges(6,3) + m(5) * t166 - mrSges(5,3) - m(4) * pkin(6) - mrSges(4,3);
t154 = qJD(3) + qJD(4);
t286 = -t277 * t111 + t276 * t154 + t289;
t285 = t280 * t112 + t278 * t154 - t290;
t196 = qJD(3) * t166;
t127 = t161 * t196;
t128 = t164 * t196;
t165 = cos(qJ(2));
t178 = t121 * t165;
t137 = t166 * t161;
t138 = t166 * t164;
t73 = t160 * t137 - t163 * t138;
t272 = qJD(1) * t178 - qJD(4) * t73 - t127 * t160 + t163 * t128;
t206 = qJD(4) * t163;
t207 = qJD(4) * t160;
t212 = qJD(1) * t165;
t271 = t163 * t127 + t160 * t128 + t137 * t206 + t138 * t207 + t183 * t212;
t158 = sin(pkin(8));
t159 = cos(pkin(8));
t284 = g(1) * t159 + g(2) * t158;
t202 = qJD(2) * qJD(3);
t129 = qJDD(2) * t164 - t161 * t202;
t283 = t129 / 0.2e1;
t176 = t183 * qJD(4);
t67 = -qJD(3) * t183 - t176;
t282 = -qJ(5) * t67 - qJD(5) * t121 + t272;
t177 = t121 * qJD(4);
t68 = -qJD(3) * t121 - t177;
t281 = qJ(5) * t68 - qJD(5) * t183 + t271;
t152 = qJDD(3) + qJDD(4);
t130 = qJDD(2) * t161 + t164 * t202;
t45 = -qJD(2) * t177 + t129 * t163 - t130 * t160;
t29 = -mrSges(6,2) * t152 + mrSges(6,3) * t45;
t30 = -mrSges(5,2) * t152 + mrSges(5,3) * t45;
t275 = t29 + t30;
t231 = mrSges(6,3) * t111;
t76 = -mrSges(6,2) * t154 - t231;
t233 = mrSges(5,3) * t111;
t77 = -mrSges(5,2) * t154 - t233;
t274 = t77 + t76;
t228 = t112 * mrSges(6,3);
t78 = mrSges(6,1) * t154 - t228;
t232 = mrSges(5,3) * t112;
t79 = mrSges(5,1) * t154 - t232;
t273 = -t79 - t78;
t261 = m(5) * pkin(3);
t270 = -mrSges(4,1) - t261;
t106 = t112 * qJ(5);
t162 = sin(qJ(2));
t205 = t162 * qJD(1);
t139 = qJD(2) * pkin(6) + t205;
t191 = pkin(7) * qJD(2) + t139;
t101 = t191 * t164;
t88 = t160 * t101;
t100 = t191 * t161;
t91 = qJD(3) * pkin(3) - t100;
t46 = t163 * t91 - t88;
t14 = -t106 + t46;
t103 = t183 * t162;
t211 = qJD(2) * t161;
t135 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t211;
t204 = t164 * qJD(2);
t136 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t204;
t269 = t164 * t135 + t161 * t136;
t268 = t164 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t129) - t161 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t130);
t267 = -t135 * t161 + t136 * t164;
t203 = qJD(1) * qJD(2);
t142 = t165 * t203;
t132 = t162 * qJDD(1) + t142;
t117 = qJDD(2) * pkin(6) + t132;
t209 = qJD(3) * t161;
t65 = t164 * t117 - t139 * t209;
t208 = qJD(3) * t164;
t66 = -t117 * t161 - t139 * t208;
t184 = -t161 * t66 + t164 * t65;
t266 = mrSges(5,1) * t150 + t291 * t151;
t219 = t159 * t165;
t97 = -t150 * t219 + t151 * t158;
t265 = -t291 * (-t150 * t158 - t151 * t219) - t293 * t97;
t220 = t158 * t165;
t95 = -t150 * t220 - t151 * t159;
t264 = -t291 * (t150 * t159 - t151 * t220) - t293 * t95;
t62 = mrSges(6,1) * t111 + mrSges(6,2) * t112;
t63 = mrSges(5,1) * t111 + mrSges(5,2) * t112;
t263 = t189 * qJD(2) + t62 + t63;
t167 = qJD(2) ^ 2;
t260 = m(6) * pkin(4);
t252 = t112 / 0.2e1;
t245 = pkin(3) * t161;
t244 = pkin(3) * t163;
t242 = pkin(4) * t112;
t239 = g(3) * t162;
t50 = -t163 * t100 - t88;
t230 = Ifges(4,4) * t161;
t229 = Ifges(4,4) * t164;
t223 = qJ(5) * t111;
t216 = t161 * t165;
t90 = t163 * t101;
t213 = t164 * t165;
t210 = qJD(2) * t162;
t201 = pkin(3) * t211;
t200 = pkin(3) * t209;
t44 = -qJD(2) * t176 + t129 * t160 + t130 * t163;
t9 = -t45 * mrSges(6,1) + t44 * mrSges(6,2);
t141 = t162 * t203;
t49 = t100 * t160 - t90;
t72 = t163 * t137 + t138 * t160;
t131 = qJDD(1) * t165 - t141;
t188 = mrSges(4,1) * t161 + mrSges(4,2) * t164;
t186 = t164 * Ifges(4,2) + t230;
t185 = Ifges(4,5) * t164 - Ifges(4,6) * t161;
t47 = t160 * t91 + t90;
t43 = qJDD(3) * pkin(3) - pkin(7) * t130 + t66;
t48 = pkin(7) * t129 + t65;
t5 = -t101 * t207 + t160 * t43 + t163 * t48 + t91 * t206;
t140 = -qJD(2) * pkin(2) - t212;
t180 = t140 * t188;
t179 = t161 * (Ifges(4,1) * t164 - t230);
t116 = -qJDD(2) * pkin(2) - t131;
t172 = t140 * t162 + (t161 ^ 2 + t164 ^ 2) * t165 * t139;
t118 = -qJD(2) * t145 - t212;
t69 = -pkin(3) * t129 + t116;
t6 = -qJD(4) * t47 - t160 * t48 + t163 * t43;
t13 = pkin(4) * t154 + t14;
t15 = t47 - t223;
t2 = pkin(4) * t152 - qJ(5) * t44 - qJD(5) * t112 + t6;
t3 = qJ(5) * t45 - qJD(5) * t111 + t5;
t64 = pkin(4) * t111 + qJD(5) + t118;
t168 = -t5 * mrSges(5,2) + t6 * mrSges(5,1) + t2 * mrSges(6,1) - t3 * mrSges(6,2) - t64 * (mrSges(6,1) * t112 - mrSges(6,2) * t111) - t118 * (mrSges(5,1) * t112 - mrSges(5,2) * t111) + t15 * t228 - t13 * t231 + t47 * t232 - t46 * t233 + t276 * t45 + t278 * t44 - (-t280 * t111 - t289) * t112 / 0.2e1 + t286 * t252 - (-t111 * t278 - t112 * t276) * t154 / 0.2e1 + (Ifges(5,3) + Ifges(6,3)) * t152 + (-t277 * t112 + t285 - t290) * t111 / 0.2e1;
t146 = Ifges(4,4) * t204;
t144 = pkin(4) + t244;
t133 = -pkin(4) * t150 - t245;
t110 = Ifges(4,1) * t211 + Ifges(4,5) * qJD(3) + t146;
t109 = Ifges(4,6) * qJD(3) + qJD(2) * t186;
t102 = t121 * t162;
t94 = pkin(4) * t183 - t145;
t75 = t201 + t242;
t70 = -mrSges(4,1) * t129 + mrSges(4,2) * t130;
t59 = -pkin(4) * t68 + t200;
t58 = -qJ(5) * t183 + t73;
t57 = -qJ(5) * t121 + t72;
t28 = mrSges(5,1) * t152 - mrSges(5,3) * t44;
t27 = mrSges(6,1) * t152 - mrSges(6,3) * t44;
t20 = -qJD(2) * t178 + t103 * t154;
t19 = -t111 * t165 + t162 * t68;
t18 = -t106 + t50;
t17 = t49 + t223;
t12 = -pkin(4) * t45 + qJDD(5) + t69;
t10 = -mrSges(5,1) * t45 + mrSges(5,2) * t44;
t1 = [m(2) * qJDD(1) - t273 * t20 + t274 * t19 - t275 * t103 - (t27 + t28) * t102 + (-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3) + (qJDD(2) * mrSges(3,1) - t167 * mrSges(3,2) + qJD(2) * t267 - t10 - t70 - t9) * t165 + (-t167 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + qJD(2) * t263 - qJD(3) * t269 + t268) * t162 + m(4) * (qJD(2) * t172 - t116 * t165 + t162 * t184) + m(3) * (t131 * t165 + t132 * t162) + m(5) * (-t102 * t6 - t103 * t5 + t118 * t210 - t165 * t69 + t19 * t47 + t20 * t46) + m(6) * (-t102 * t2 - t103 * t3 - t12 * t165 + t13 * t20 + t15 * t19 + t210 * t64); t184 * mrSges(4,3) - t267 * t212 + (m(4) * t184 - t135 * t208 - t136 * t209 + t268) * pkin(6) + t285 * t67 / 0.2e1 + t286 * t68 / 0.2e1 + t110 * t208 / 0.2e1 - t109 * t209 / 0.2e1 + t116 * t189 + (t180 + t185 * qJD(3) / 0.2e1) * qJD(3) + t130 * t229 / 0.2e1 + (t131 + t141) * mrSges(3,1) + qJDD(3) * (Ifges(4,5) * t161 + Ifges(4,6) * t164) + (-t13 * t67 + t15 * t68) * mrSges(6,3) + (-t46 * t67 + t47 * t68) * mrSges(5,3) - t145 * t10 + t118 * (-mrSges(5,1) * t68 + mrSges(5,2) * t67) + t94 * t9 + t64 * (-mrSges(6,1) * t68 + mrSges(6,2) * t67) - pkin(2) * t70 + t72 * t28 + t73 * t30 + t59 * t62 + t57 * t27 + t58 * t29 + (-m(5) * t118 - m(6) * t64 - t263) * t205 + (-pkin(2) * t116 - qJD(1) * t172) * m(4) + (t287 * g(3) + t284 * t288) * t162 + (-t288 * g(3) + t284 * t287) * t165 + (t164 * (-Ifges(4,2) * t161 + t229) + t179) * t202 / 0.2e1 + t164 * (Ifges(4,4) * t130 + Ifges(4,2) * t129) / 0.2e1 + t281 * t76 + t282 * t78 + (t12 * t94 + t13 * t282 + t15 * t281 + t2 * t57 + t3 * t58 + t59 * t64) * m(6) + (Ifges(4,1) * t130 + Ifges(4,4) * t283) * t161 - (t277 * t68 + t279 * t67) * t111 / 0.2e1 + (t279 * t68 + t280 * t67) * t252 + t271 * t77 + t272 * t79 + (t118 * t200 - t145 * t69 + t271 * t47 + t272 * t46 + t5 * t73 + t6 * t72) * m(5) + (t276 * t68 + t278 * t67) * t154 / 0.2e1 + (-t132 + t142) * mrSges(3,2) + t186 * t283 + t63 * t200 + (mrSges(5,1) * t69 + mrSges(6,1) * t12 - mrSges(5,3) * t5 - mrSges(6,3) * t3 - t152 * t276 - t277 * t45 - t279 * t44) * t183 + Ifges(3,3) * qJDD(2) + (mrSges(5,2) * t69 + mrSges(6,2) * t12 - mrSges(5,3) * t6 - mrSges(6,3) * t2 + t152 * t278 + t279 * t45 + t280 * t44) * t121; (m(5) * t245 - m(6) * t133 + mrSges(6,1) * t150 + t188 + t266) * t239 + t269 * t139 + t109 * t211 / 0.2e1 - t185 * t202 / 0.2e1 - t63 * t201 - m(5) * (t118 * t201 + t46 * t49 + t47 * t50) + t144 * t27 + Ifges(4,6) * t129 + Ifges(4,5) * t130 - t17 * t78 - t49 * t79 - t75 * t62 - t18 * t76 - t50 * t77 - t65 * mrSges(4,2) + t66 * mrSges(4,1) + t168 - qJD(2) * t180 - t167 * t179 / 0.2e1 - (-Ifges(4,2) * t211 + t110 + t146) * t204 / 0.2e1 + (-(-t158 * t213 + t159 * t161) * mrSges(4,2) - m(6) * (t133 * t220 - t134 * t159) + t270 * (-t158 * t216 - t159 * t164) + t264) * g(2) + (-(-t158 * t161 - t159 * t213) * mrSges(4,2) - m(6) * (t133 * t219 + t134 * t158) + t270 * (t158 * t164 - t159 * t216) + t265) * g(1) + (t160 * t5 + t163 * t6 + (-t160 * t46 + t163 * t47) * qJD(4)) * t261 + t28 * t244 + (-t13 * t17 + t144 * t2 - t15 * t18 - t64 * t75) * m(6) + ((t160 * t3 + (-t13 * t160 + t15 * t163) * qJD(4)) * m(6) + t273 * t207 + t274 * t206 + t275 * t160) * pkin(3) + Ifges(4,3) * qJDD(3); t15 * t78 + t47 * t79 - t14 * t76 - t46 * t77 + pkin(4) * t27 + t168 - t62 * t242 + t2 * t260 - m(6) * (t64 * t242 + (-t13 + t14) * t15) + (-(-mrSges(6,1) - t260) * t150 + t266) * t239 + (-t260 * t95 + t264) * g(2) + (-t260 * t97 + t265) * g(1); t111 * t76 + t112 * t78 + (g(3) * t165 + t15 * t111 + t13 * t112 - t284 * t162 + t12) * m(6) + t9;];
tau = t1;
