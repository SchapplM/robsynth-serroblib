% Calculate vector of inverse dynamics joint torques for
% S5RRPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:54
% EndTime: 2019-12-31 19:29:12
% DurationCPUTime: 9.88s
% Computational Cost: add. (3375->482), mult. (7686->613), div. (0->0), fcn. (5153->10), ass. (0->223)
t162 = -qJ(3) - pkin(6);
t164 = sin(qJ(2));
t134 = t162 * t164;
t126 = qJD(1) * t134;
t167 = cos(qJ(2));
t136 = t162 * t167;
t127 = qJD(1) * t136;
t161 = sin(pkin(8));
t225 = t161 * t127;
t232 = cos(pkin(8));
t78 = t232 * t126 + t225;
t295 = -t78 + qJD(4);
t294 = mrSges(4,1) + mrSges(5,1);
t288 = Ifges(4,1) + Ifges(5,1);
t286 = Ifges(5,4) + Ifges(4,5);
t122 = t161 * t167 + t164 * t232;
t108 = t122 * qJD(1);
t244 = t108 * pkin(7);
t293 = -t244 + t295;
t165 = sin(qJ(1));
t292 = g(2) * t165;
t287 = -Ifges(4,4) + Ifges(5,5);
t285 = Ifges(5,6) - Ifges(4,6);
t154 = t167 * pkin(2);
t284 = t154 + pkin(1);
t131 = -qJD(1) * t284 + qJD(3);
t291 = -qJ(4) * t108 + t131;
t135 = -mrSges(3,1) * t167 + mrSges(3,2) * t164;
t157 = qJ(2) + pkin(8);
t152 = sin(t157);
t153 = cos(t157);
t290 = t135 - t294 * t153 - (-mrSges(4,2) + mrSges(5,3)) * t152;
t218 = qJD(1) * qJD(2);
t206 = t164 * t218;
t217 = qJDD(1) * t167;
t128 = -t206 + t217;
t129 = qJDD(1) * t164 + t167 * t218;
t80 = t161 * t128 + t129 * t232;
t263 = t80 / 0.2e1;
t289 = t128 / 0.2e1;
t242 = qJD(2) / 0.2e1;
t202 = t232 * t167;
t224 = qJD(1) * t164;
t106 = -qJD(1) * t202 + t161 * t224;
t101 = Ifges(4,4) * t106;
t234 = Ifges(5,5) * t106;
t283 = t286 * qJD(2) + t288 * t108 - t101 + t234;
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t247 = pkin(7) * t106;
t203 = t232 * t127;
t77 = t126 * t161 - t203;
t37 = t77 + t247;
t209 = t232 * pkin(2);
t148 = -t209 - pkin(3);
t143 = -pkin(4) + t148;
t250 = pkin(2) * t161;
t145 = qJ(4) + t250;
t93 = t143 * t166 - t145 * t163;
t282 = qJD(5) * t93 - t163 * t37 + t166 * t293;
t94 = t143 * t163 + t145 * t166;
t281 = -qJD(5) * t94 - t163 * t293 - t166 * t37;
t239 = mrSges(4,3) * t106;
t241 = mrSges(5,2) * t106;
t92 = qJD(2) * mrSges(5,3) - t241;
t280 = -qJD(2) * mrSges(4,2) - t239 + t92;
t238 = mrSges(4,3) * t108;
t240 = mrSges(5,2) * t108;
t279 = qJD(2) * t294 - t238 - t240;
t233 = qJDD(2) / 0.2e1;
t168 = cos(qJ(1));
t273 = g(1) * t168 + t292;
t278 = t273 * t152;
t226 = t153 * t168;
t227 = t152 * t168;
t277 = pkin(3) * t226 + qJ(4) * t227;
t249 = pkin(2) * t164;
t261 = -pkin(3) - pkin(4);
t276 = -m(5) * (-pkin(3) * t152 - t249) + t152 * mrSges(5,1) - t153 * mrSges(5,3) - m(6) * (t152 * t261 - t249);
t150 = pkin(6) * t217;
t119 = -pkin(6) * t206 + t150;
t120 = t129 * pkin(6);
t274 = t119 * t167 + t120 * t164;
t156 = -qJD(2) + qJD(5);
t271 = 0.2e1 * t233;
t270 = -m(3) * pkin(1) - mrSges(2,1) + t290;
t269 = mrSges(2,2) - mrSges(5,2) - mrSges(4,3) - m(3) * pkin(6) - mrSges(3,3) - m(6) * (-pkin(7) - t162) + mrSges(6,3);
t268 = (-g(1) * t226 - t153 * t292) * qJ(4);
t55 = t106 * t166 - t108 * t163;
t266 = -t55 / 0.2e1;
t184 = t106 * t163 + t108 * t166;
t265 = -t184 / 0.2e1;
t264 = t184 / 0.2e1;
t262 = m(5) + m(6);
t118 = qJD(2) * pkin(2) + t126;
t69 = t118 * t232 + t225;
t179 = qJD(4) - t69;
t30 = qJD(2) * t261 + t179 - t244;
t70 = t161 * t118 - t203;
t58 = qJD(2) * qJ(4) + t70;
t35 = t58 + t247;
t6 = -t163 * t35 + t166 * t30;
t260 = t6 * mrSges(6,3);
t7 = t163 * t30 + t166 * t35;
t259 = t7 * mrSges(6,3);
t258 = -t106 / 0.2e1;
t257 = t106 / 0.2e1;
t255 = t108 / 0.2e1;
t252 = -t156 / 0.2e1;
t251 = Ifges(6,4) * t184;
t248 = pkin(6) * t167;
t146 = t153 * pkin(3);
t221 = qJD(3) * t164;
t65 = qJDD(2) * pkin(2) - qJ(3) * t129 - qJD(1) * t221 - t120;
t222 = qJD(2) * t164;
t213 = pkin(6) * t222;
t220 = qJD(3) * t167;
t73 = qJ(3) * t128 + t150 + (-t213 + t220) * qJD(1);
t25 = t161 * t65 + t232 * t73;
t204 = qJD(2) * t162;
t104 = t164 * t204 + t220;
t105 = t167 * t204 - t221;
t44 = t232 * t104 + t161 * t105;
t237 = Ifges(3,4) * t164;
t236 = Ifges(3,4) * t167;
t235 = Ifges(4,4) * t108;
t230 = qJDD(1) * pkin(1);
t144 = t152 * qJ(4);
t84 = t161 * t134 - t232 * t136;
t223 = qJD(1) * t167;
t216 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t224) * t248;
t215 = pkin(2) * t224;
t214 = pkin(2) * t222;
t211 = t146 + t144 + t154;
t79 = -t128 * t232 + t129 * t161;
t208 = t79 * mrSges(4,1) + t80 * mrSges(4,2);
t207 = t79 * mrSges(5,1) - t80 * mrSges(5,3);
t43 = t104 * t161 - t232 * t105;
t61 = -qJDD(2) * mrSges(5,1) + t80 * mrSges(5,2);
t200 = -t284 - t144;
t83 = -t232 * t134 - t136 * t161;
t141 = t168 * t284;
t199 = -t165 * t162 + t141;
t21 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t25;
t198 = qJ(4) * t122 + t284;
t12 = qJD(5) * t55 + t163 * t79 + t166 * t80;
t13 = -qJD(5) * t184 - t163 * t80 + t166 * t79;
t196 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t103 = t152 * t166 - t153 * t163;
t85 = t103 * t165;
t183 = t152 * t163 + t153 * t166;
t86 = t183 * t165;
t195 = t85 * mrSges(6,1) - t86 * mrSges(6,2);
t87 = t163 * t226 - t166 * t227;
t88 = t183 * t168;
t194 = -t87 * mrSges(6,1) - t88 * mrSges(6,2);
t24 = -t161 * t73 + t232 * t65;
t191 = mrSges(3,1) * t164 + mrSges(3,2) * t167;
t188 = -mrSges(6,1) * t183 - t103 * mrSges(6,2);
t187 = t167 * Ifges(3,2) + t237;
t186 = Ifges(3,5) * t167 - Ifges(3,6) * t164;
t39 = -mrSges(6,2) * t156 + mrSges(6,3) * t55;
t40 = mrSges(6,1) * t156 - mrSges(6,3) * t184;
t185 = -t163 * t40 + t166 * t39;
t45 = -pkin(7) * t122 + t83;
t176 = -t161 * t164 + t202;
t46 = -pkin(7) * t176 + t84;
t17 = -t163 * t46 + t166 * t45;
t18 = t163 * t45 + t166 * t46;
t71 = -t122 * t163 - t166 * t176;
t72 = t122 * t166 - t163 * t176;
t182 = qJDD(4) - t24;
t181 = -qJ(4) * t106 - t215;
t180 = pkin(1) * t191;
t178 = pkin(2) * t128 - qJDD(3) + t230;
t177 = t164 * (Ifges(3,1) * t167 - t237);
t109 = t176 * qJD(2);
t172 = qJ(4) * t109 + qJD(4) * t122 - t214;
t170 = qJ(4) * t80 + qJD(4) * t108 + t178;
t10 = -t80 * pkin(7) + qJDD(2) * t261 + t182;
t11 = pkin(7) * t79 + t21;
t1 = qJD(5) * t6 + t10 * t163 + t11 * t166;
t155 = -qJDD(2) + qJDD(5);
t2 = -qJD(5) * t7 + t10 * t166 - t11 * t163;
t169 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t12 + Ifges(6,6) * t13 + Ifges(6,3) * t155;
t151 = Ifges(3,4) * t223;
t133 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t223;
t116 = Ifges(3,1) * t224 + Ifges(3,5) * qJD(2) + t151;
t115 = Ifges(3,6) * qJD(2) + qJD(1) * t187;
t107 = t122 * qJD(2);
t100 = Ifges(5,5) * t108;
t68 = -pkin(3) * t176 - t198;
t67 = mrSges(4,1) * t106 + mrSges(4,2) * t108;
t66 = mrSges(5,1) * t106 - mrSges(5,3) * t108;
t62 = -mrSges(5,2) * t79 + qJDD(2) * mrSges(5,3);
t60 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t80;
t59 = -qJDD(2) * mrSges(4,2) - mrSges(4,3) * t79;
t51 = -t106 * Ifges(4,2) + Ifges(4,6) * qJD(2) + t235;
t50 = Ifges(5,6) * qJD(2) + t106 * Ifges(5,3) + t100;
t49 = -qJD(2) * pkin(3) + t179;
t47 = Ifges(6,4) * t55;
t42 = pkin(3) * t108 - t181;
t41 = pkin(3) * t106 + t291;
t36 = -t176 * t261 + t198;
t34 = pkin(3) * t107 - t172;
t33 = pkin(7) * t107 + t44;
t32 = -pkin(7) * t109 + t43;
t31 = t108 * t261 + t181;
t29 = t106 * t261 - t291;
t28 = t107 * t261 + t172;
t27 = -qJD(5) * t72 + t107 * t166 - t109 * t163;
t26 = qJD(5) * t71 + t107 * t163 + t109 * t166;
t23 = -qJDD(2) * pkin(3) + t182;
t22 = -mrSges(6,1) * t55 + mrSges(6,2) * t184;
t20 = Ifges(6,1) * t184 + t156 * Ifges(6,5) + t47;
t19 = t55 * Ifges(6,2) + t156 * Ifges(6,6) + t251;
t16 = pkin(3) * t79 - t170;
t9 = -mrSges(6,2) * t155 + mrSges(6,3) * t13;
t8 = mrSges(6,1) * t155 - mrSges(6,3) * t12;
t5 = t261 * t79 + t170;
t4 = -qJD(5) * t18 - t163 * t33 + t166 * t32;
t3 = qJD(5) * t17 + t163 * t32 + t166 * t33;
t14 = [t167 * (Ifges(3,4) * t129 + Ifges(3,2) * t128 + Ifges(3,6) * qJDD(2)) / 0.2e1 + (-m(4) * t24 + m(5) * t23 - t60 + t61) * t83 + Ifges(3,6) * t167 * t233 + t36 * t196 + t67 * t214 + t156 * (Ifges(6,5) * t26 + Ifges(6,6) * t27) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t128 + mrSges(3,2) * t129) + (t283 / 0.2e1 + mrSges(4,2) * t131 + t49 * mrSges(5,2) - t69 * mrSges(4,3) - mrSges(5,3) * t41 + Ifges(4,4) * t258 + Ifges(5,5) * t257 + t242 * t286 + t255 * t288) * t109 + t34 * t66 + t55 * (Ifges(6,4) * t26 + Ifges(6,2) * t27) / 0.2e1 + (-t178 * mrSges(4,2) + t23 * mrSges(5,2) - t24 * mrSges(4,3) - t16 * mrSges(5,3) + (-Ifges(4,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t79 + t288 * t263 + t286 * t233) * t122 + m(4) * (t131 * t214 + t178 * t284) + (t86 * mrSges(6,1) + t85 * mrSges(6,2) + ((m(4) + m(5)) * t162 + t269) * t168 + (m(4) * t284 - m(5) * (t200 - t146) - m(6) * (t153 * t261 + t200) - t270) * t165) * g(1) - t284 * t208 + (m(4) * t25 + m(5) * t21 + t59 + t62) * t84 + (-mrSges(6,1) * t5 + mrSges(6,3) * t1 + Ifges(6,4) * t12 + Ifges(6,2) * t13 + Ifges(6,6) * t155) * t71 + (m(4) * t70 + m(5) * t58 + t280) * t44 + (Ifges(5,3) * t257 - t58 * mrSges(5,2) - t70 * mrSges(4,3) + t131 * mrSges(4,1) + t41 * mrSges(5,1) + t50 / 0.2e1 - t51 / 0.2e1 - Ifges(4,2) * t258 + t287 * t255 + t285 * t242) * t107 + (qJDD(2) * t286 + t287 * t79 + t288 * t80) * t122 / 0.2e1 + t129 * t236 / 0.2e1 + t3 * t39 + t4 * t40 + t28 * t22 + t29 * (-mrSges(6,1) * t27 + mrSges(6,2) * t26) + t26 * t20 / 0.2e1 + t27 * t19 / 0.2e1 + t17 * t8 + t18 * t9 + t27 * t259 + (Ifges(6,1) * t26 + Ifges(6,4) * t27) * t264 + (t128 * t248 + t274) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t274) + (-m(5) * (t199 + t277) - m(4) * t199 - m(6) * (pkin(4) * t226 + t141 + t277) - t88 * mrSges(6,1) + t87 * mrSges(6,2) + t270 * t168 + t269 * t165) * g(2) + (-m(4) * t69 + m(5) * t49 - t279) * t43 + (t186 * t242 - t216) * qJD(2) + (mrSges(6,2) * t5 - mrSges(6,3) * t2 + Ifges(6,1) * t12 + Ifges(6,4) * t13 + Ifges(6,5) * t155) * t72 + (t167 * (-Ifges(3,2) * t164 + t236) + t177) * t218 / 0.2e1 + (-pkin(6) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t129) + Ifges(3,1) * t129 + Ifges(3,4) * t289 + t271 * Ifges(3,5)) * t164 + t187 * t289 + m(6) * (t1 * t18 + t17 * t2 + t28 * t29 + t3 * t7 + t36 * t5 + t4 * t6) + t167 * t116 * t242 + t68 * t207 - t133 * t213 - t180 * t218 - t115 * t222 / 0.2e1 - t135 * t230 - t26 * t260 - qJDD(2) * mrSges(3,2) * t248 - (-mrSges(4,1) * t178 + mrSges(5,1) * t16 - mrSges(5,2) * t21 - mrSges(4,3) * t25 + (Ifges(4,2) + Ifges(5,3)) * t79 + 0.2e1 * t287 * t263 + t285 * t271) * t176 + Ifges(2,3) * qJDD(1) + m(5) * (t16 * t68 + t34 * t41); ((t161 * t25 + t232 * t24) * pkin(2) - t131 * t215 + t69 * t77 - t70 * t78) * m(4) + t59 * t250 + t51 * t255 + t70 * t238 + t58 * t240 + t49 * t241 + t145 * t62 + t148 * t61 + Ifges(3,5) * t129 - t131 * (mrSges(4,1) * t108 - mrSges(4,2) * t106) + Ifges(3,6) * t128 - t119 * mrSges(3,2) - t120 * mrSges(3,1) - t41 * (mrSges(5,1) * t108 + mrSges(5,3) * t106) + qJD(4) * t92 + t93 * t8 + t94 * t9 - t42 * t66 - (Ifges(6,5) * t252 - t29 * mrSges(6,2) - t20 / 0.2e1 + t260 + Ifges(6,1) * t265 + Ifges(6,4) * t266) * t55 - t280 * t78 + t281 * t40 + (t1 * t94 + t2 * t93 + t281 * t6 + t282 * t7 - t29 * t31 + t268) * m(6) + t282 * t39 + (-Ifges(4,2) * t108 - t101 + t283) * t257 + t285 * t79 + t286 * t80 - (-t106 * t286 + t108 * t285) * qJD(2) / 0.2e1 - (-t106 * t288 + t100 - t235 + t50) * t108 / 0.2e1 + (m(4) * t249 + mrSges(4,1) * t152 + mrSges(4,2) * t153 + t191) * t273 - t31 * t22 - t23 * mrSges(5,1) + t24 * mrSges(4,1) - t25 * mrSges(4,2) + t21 * mrSges(5,3) + (Ifges(5,3) * t108 - t234) * t258 + (t216 + (t180 - t177 / 0.2e1) * qJD(1)) * qJD(1) + (t168 * t276 + t194) * g(1) + (t165 * t276 + t195) * g(2) + t279 * t77 - (-Ifges(3,2) * t224 + t116 + t151) * t223 / 0.2e1 + (Ifges(6,6) * t252 - t19 / 0.2e1 + t29 * mrSges(6,1) + Ifges(6,4) * t265 + Ifges(6,2) * t266 - t259) * t184 + (Ifges(5,2) + Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (pkin(6) * t133 + t115 / 0.2e1) * t224 + t60 * t209 - t169 + (t145 * t21 + t148 * t23 + t295 * t58 - t41 * t42 - t49 * t77 + t268) * m(5) - t67 * t215 - t186 * t218 / 0.2e1 - t69 * t239 + (-m(6) * (pkin(4) * t153 + t211) + t188 - m(5) * t211 - m(4) * t154 + t290) * g(3); t55 * t39 - t184 * t40 + t279 * t108 + t280 * t106 - t196 + t207 + t208 + (-g(1) * t165 + g(2) * t168) * (m(4) + t262) + (-t184 * t6 + t55 * t7 - t5) * m(6) + (t106 * t58 - t108 * t49 + t16) * m(5) + (t106 * t70 + t108 * t69 - t178) * m(4); t163 * t9 + t166 * t8 + (-t22 + t66) * t108 + t185 * qJD(5) + t262 * t153 * g(3) + (-t185 - t92) * qJD(2) + t61 + (t1 * t163 - t108 * t29 - t278 + t166 * t2 + t156 * (-t163 * t6 + t166 * t7)) * m(6) + (-qJD(2) * t58 + t108 * t41 + t23 - t278) * m(5); -t29 * (mrSges(6,1) * t184 + mrSges(6,2) * t55) + (Ifges(6,1) * t55 - t251) * t265 + t19 * t264 + (Ifges(6,5) * t55 - Ifges(6,6) * t184) * t252 - t6 * t39 + t7 * t40 - g(1) * t194 - g(2) * t195 - g(3) * t188 + (t184 * t7 + t55 * t6) * mrSges(6,3) + t169 + (-Ifges(6,2) * t184 + t20 + t47) * t266;];
tau = t14;
