% Calculate time derivative of joint inertia matrix for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:37
% EndTime: 2019-03-09 16:43:47
% DurationCPUTime: 4.39s
% Computational Cost: add. (4013->347), mult. (8791->475), div. (0->0), fcn. (7656->6), ass. (0->160)
t287 = Ifges(6,1) + Ifges(7,1);
t279 = Ifges(7,4) + Ifges(6,5);
t146 = sin(qJ(5));
t149 = cos(qJ(5));
t231 = Ifges(7,5) * t149;
t233 = Ifges(6,4) * t149;
t286 = t287 * t146 - t231 + t233;
t285 = -Ifges(7,2) - Ifges(6,3);
t284 = -Ifges(6,6) + Ifges(7,6);
t147 = sin(qJ(3));
t214 = qJD(3) * t147;
t204 = pkin(2) * t214;
t216 = t146 ^ 2 + t149 ^ 2;
t283 = t204 * t216;
t232 = Ifges(7,5) * t146;
t234 = Ifges(6,4) * t146;
t282 = t287 * t149 + t232 - t234;
t281 = 2 * qJD(4);
t280 = mrSges(5,2) - mrSges(4,1);
t148 = sin(qJ(2));
t243 = cos(qJ(3));
t244 = cos(qJ(2));
t110 = t147 * t148 - t243 * t244;
t213 = qJD(5) * t146;
t111 = t147 * t244 + t148 * t243;
t258 = qJD(2) + qJD(3);
t79 = t258 * t111;
t228 = t149 * t79;
t164 = t110 * t213 - t228;
t212 = qJD(5) * t149;
t229 = t146 * t79;
t165 = t110 * t212 + t229;
t78 = t110 * t258;
t278 = -t279 * t78 + t287 * t165 + (-Ifges(6,4) + Ifges(7,5)) * t164;
t277 = t286 * t110 + t279 * t111;
t276 = (mrSges(7,2) + mrSges(6,3)) * t216;
t275 = t286 * qJD(5);
t205 = t243 * pkin(2);
t138 = -t205 - pkin(3);
t133 = -pkin(9) + t138;
t273 = -t133 * t213 + t149 * t204;
t139 = -pkin(2) * t244 - pkin(1);
t162 = -t111 * qJ(4) + t139;
t246 = pkin(3) + pkin(9);
t50 = t110 * t246 + t162;
t269 = -pkin(8) - pkin(7);
t128 = t269 * t148;
t129 = t269 * t244;
t90 = -t243 * t128 - t129 * t147;
t59 = pkin(4) * t111 + t90;
t236 = t146 * t59 + t149 * t50;
t24 = -t146 * t50 + t149 * t59;
t215 = qJD(2) * t148;
t142 = pkin(2) * t215;
t163 = qJ(4) * t78 - qJD(4) * t111 + t142;
t17 = t246 * t79 + t163;
t259 = qJD(5) * t236;
t159 = qJD(2) * t129;
t183 = qJD(2) * t128;
t91 = t147 * t128 - t129 * t243;
t39 = qJD(3) * t91 + t147 * t183 - t243 * t159;
t28 = -t78 * pkin(4) + t39;
t5 = -t146 * t17 + t149 * t28 - t259;
t3 = pkin(5) * t78 - t5;
t240 = t149 * t3;
t18 = qJ(6) * t111 + t236;
t19 = -pkin(5) * t111 - t24;
t4 = t146 * t28 + t149 * t17 + t59 * t212 - t213 * t50;
t1 = -qJ(6) * t78 + qJD(6) * t111 + t4;
t241 = t1 * t146;
t255 = -qJD(5) * (t146 * t19 + t149 * t18) - t241;
t239 = t4 * t146;
t260 = t5 * t149 + t239;
t30 = mrSges(6,2) * t78 - mrSges(6,3) * t164;
t31 = -mrSges(6,1) * t78 - mrSges(6,3) * t165;
t32 = t78 * mrSges(7,1) + mrSges(7,2) * t165;
t33 = -mrSges(7,2) * t164 - mrSges(7,3) * t78;
t272 = m(6) * ((-t146 * t24 + t149 * t236) * qJD(5) + t260) + m(7) * (-t240 - t255) + (t31 - t32) * t149 + (t30 + t33) * t146;
t34 = pkin(3) * t79 + t163;
t271 = -0.2e1 * t34;
t270 = 0.2e1 * t142;
t267 = Ifges(3,1) - Ifges(3,2);
t226 = t110 * t146;
t68 = mrSges(6,1) * t111 - mrSges(6,3) * t226;
t69 = -mrSges(7,1) * t111 + mrSges(7,2) * t226;
t264 = t68 - t69;
t225 = t110 * t149;
t70 = -mrSges(6,2) * t111 + mrSges(6,3) * t225;
t71 = mrSges(7,2) * t225 + mrSges(7,3) * t111;
t263 = t70 + t71;
t123 = t146 * mrSges(6,1) + t149 * mrSges(6,2);
t262 = mrSges(5,3) + t123;
t261 = t284 * t164 + t279 * t165 + t285 * t78;
t171 = t146 * pkin(5) - qJ(6) * t149;
t252 = 2 * m(5);
t251 = 0.2e1 * m(6);
t250 = 0.2e1 * m(7);
t179 = t149 * mrSges(6,1) - t146 * mrSges(6,2);
t114 = t179 * qJD(5);
t249 = 0.2e1 * t114;
t122 = t146 * mrSges(7,1) - t149 * mrSges(7,3);
t248 = 0.2e1 * t122;
t247 = -t78 / 0.2e1;
t242 = pkin(2) * t147;
t237 = t78 * mrSges(5,1);
t230 = pkin(2) * qJD(3);
t193 = qJD(3) * t243;
t185 = pkin(2) * t193;
t130 = t185 + qJD(4);
t135 = qJ(4) + t242;
t224 = t130 * t135;
t125 = -Ifges(6,2) * t146 + t233;
t220 = t149 * t125;
t219 = t133 * t283;
t218 = t246 * t283;
t140 = Ifges(7,6) * t212;
t201 = mrSges(7,2) * t213;
t217 = pkin(5) * t201 + t140;
t211 = qJD(6) * t146;
t202 = t243 * mrSges(4,2);
t194 = qJD(2) * t244;
t191 = t212 / 0.2e1;
t190 = t111 * t204;
t38 = -t128 * t193 - t129 * t214 - t147 * t159 - t243 * t183;
t180 = -t38 * t91 + t39 * t90;
t119 = qJ(4) + t171;
t178 = t149 * mrSges(7,1) + t146 * mrSges(7,3);
t175 = Ifges(6,2) * t149 + t234;
t174 = Ifges(6,5) * t146 + Ifges(6,6) * t149;
t173 = -Ifges(7,3) * t149 + t232;
t172 = pkin(5) * t149 + qJ(6) * t146;
t167 = qJ(4) * t130 + qJD(4) * t135;
t166 = -pkin(4) - t172;
t93 = pkin(5) * t212 + qJ(6) * t213 - t149 * qJD(6) + qJD(4);
t158 = (-mrSges(7,2) * qJ(6) - Ifges(6,6)) * t149 - t279 * t146;
t115 = t173 * qJD(5);
t116 = t175 * qJD(5);
t124 = Ifges(7,3) * t146 + t231;
t155 = -t220 * qJD(5) + t124 * t212 - t275 * t149 + (-t282 * qJD(5) - t115 + t116) * t146;
t154 = -m(7) * t171 - t122 - t123;
t153 = m(7) * t211 + qJD(5) * t154;
t113 = t178 * qJD(5);
t12 = Ifges(7,5) * t165 - t78 * Ifges(7,6) + Ifges(7,3) * t164;
t13 = Ifges(6,4) * t165 - Ifges(6,2) * t164 - t78 * Ifges(6,6);
t27 = -pkin(4) * t79 - t38;
t35 = t110 * t166 + t91;
t46 = Ifges(7,6) * t111 + t110 * t173;
t47 = Ifges(6,6) * t111 + t110 * t175;
t60 = -t110 * pkin(4) + t91;
t7 = t166 * t79 + (qJD(5) * t171 - t211) * t110 - t38;
t152 = -t111 * t174 * qJD(5) / 0.2e1 + mrSges(7,2) * t240 + t60 * t114 + t7 * t122 + t27 * t123 + t35 * t113 - t47 * t212 / 0.2e1 + t111 * (-Ifges(7,4) * t213 + t140) / 0.2e1 - t124 * t228 / 0.2e1 + t46 * t191 + t24 * mrSges(6,3) * t213 + (Ifges(5,4) - Ifges(4,5)) * t78 + t280 * t39 + (-t116 / 0.2e1 + t115 / 0.2e1) * t225 + (-t13 / 0.2e1 + t12 / 0.2e1 + t284 * t247) * t146 + (Ifges(5,5) - Ifges(4,6) + t220 / 0.2e1) * t79 + (qJD(5) * t124 - t275) * t226 / 0.2e1 - (t110 * t125 + t277) * t213 / 0.2e1 + t282 * (t110 * t191 + t229 / 0.2e1) + (t278 / 0.2e1 + t279 * t247) * t149;
t151 = (-t146 * t264 + t149 * t263) * qJD(5) + t272;
t105 = t119 + t242;
t92 = t185 + t93;
t63 = t110 * pkin(3) + t162;
t62 = t179 * t110;
t61 = t178 * t110;
t21 = mrSges(6,1) * t164 + mrSges(6,2) * t165;
t20 = mrSges(7,1) * t164 - mrSges(7,3) * t165;
t2 = [0.2e1 * m(4) * (t139 * t142 + t180) + (-t12 + t13) * t225 + (t34 * t63 + t180) * t252 + (t1 * t18 + t19 * t3 + t35 * t7) * t250 + t165 * t277 + t278 * t226 + (t236 * t4 + t24 * t5 + t27 * t60) * t251 + 0.2e1 * t236 * t30 + t164 * (-t47 + t46) + (0.2e1 * Ifges(3,4) * t244 + t148 * t267) * t194 + (-0.2e1 * Ifges(3,4) * t148 + t244 * t267) * t215 + 0.2e1 * (Ifges(5,6) + Ifges(4,4)) * (t110 * t78 - t111 * t79) + (mrSges(4,1) * t270 + mrSges(5,2) * t271 + 0.2e1 * (Ifges(5,3) + Ifges(4,2)) * t79 + (-Ifges(7,4) * t146 + Ifges(7,6) * t149 - t174) * t78) * t110 + (mrSges(4,2) * t270 + mrSges(5,3) * t271 + (-(2 * Ifges(4,1)) - (2 * Ifges(5,2)) + t285) * t78 + t261) * t111 + 0.2e1 * (mrSges(4,3) + mrSges(5,1)) * (t110 * t38 + t111 * t39 - t78 * t90 - t79 * t91) + 0.2e1 * t139 * (mrSges(4,1) * t79 - mrSges(4,2) * t78) + 0.2e1 * t63 * (-mrSges(5,2) * t79 + mrSges(5,3) * t78) - 0.2e1 * t7 * t61 - 0.2e1 * t27 * t62 + 0.2e1 * t5 * t68 + 0.2e1 * t3 * t69 + 0.2e1 * t4 * t70 + 0.2e1 * t1 * t71 + 0.2e1 * t60 * t21 + 0.2e1 * t35 * t20 + 0.2e1 * t24 * t31 + 0.2e1 * t19 * t32 + 0.2e1 * t18 * t33 - 0.2e1 * pkin(1) * (mrSges(3,1) * t148 + mrSges(3,2) * t244) * qJD(2); t263 * t146 * t204 + m(7) * (t105 * t7 + t35 * t92 + (t146 * t18 - t149 * t19) * t204) + m(5) * (t130 * t91 - t135 * t38 + t138 * t39 + t204 * t90) - t19 * t201 + t152 + Ifges(3,5) * t194 + (-t212 * t236 - t260) * mrSges(6,3) - t130 * t62 + t135 * t21 + t105 * t20 - t92 * t61 - t38 * mrSges(5,3) + t38 * mrSges(4,2) + (-mrSges(3,1) * t194 + mrSges(3,2) * t215) * pkin(7) + (-t110 * t185 + t205 * t78 - t242 * t79 + t190) * mrSges(4,3) + (-t18 * t212 - t241) * mrSges(7,2) + (-t110 * t130 - t135 * t79 + t190) * mrSges(5,1) + m(6) * (t130 * t60 + t135 * t27 + (t146 * t236 + t149 * t24) * t204) + m(4) * (-t243 * t39 - t147 * t38 + (t147 * t90 + t243 * t91) * qJD(3)) * pkin(2) - Ifges(3,6) * t215 - t138 * t237 + t273 * t264 + (t212 * t263 + t272) * t133; 0.2e1 * t105 * t113 + t135 * t249 + t92 * t248 + 0.2e1 * t262 * t130 + t224 * t252 + (t105 * t92 + t219) * t250 + (t219 + t224) * t251 + (-0.2e1 * t202 + (t138 * t252 - 0.2e1 * mrSges(4,1) + 0.2e1 * mrSges(5,2) - 0.2e1 * t276) * t147) * t230 + t155; t255 * mrSges(7,2) + (-t239 + (-t5 - t259) * t149) * mrSges(6,3) + (-mrSges(5,3) + mrSges(4,2)) * t38 + m(6) * (qJ(4) * t27 + qJD(4) * t60) + m(7) * (t119 * t7 + t93 * t35) + t152 + m(5) * (-pkin(3) * t39 - qJ(4) * t38 + qJD(4) * t91) + (pkin(3) * t78 - qJ(4) * t79 - qJD(4) * t110) * mrSges(5,1) - t151 * t246 + t119 * t20 - t93 * t61 - qJD(4) * t62 + qJ(4) * t21; (t93 + t92) * t122 + (t135 + qJ(4)) * t114 + (t105 + t119) * t113 + m(7) * (t105 * t93 + t119 * t92 - t218) + m(6) * (t167 - t218) + m(5) * t167 + (-t202 + (-m(5) * pkin(3) - t276 + t280) * t147) * t230 + t155 + t262 * (qJD(4) + t130); t93 * t248 + t262 * t281 + t155 + 0.2e1 * (m(7) * t93 + t113) * t119 + (t249 + (m(5) + m(6)) * t281) * qJ(4); m(5) * t39 + t151 - t237; (m(5) + 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t216) * t204; 0; 0; m(7) * (-pkin(5) * t3 + qJ(6) * t1 + qJD(6) * t18) - t4 * mrSges(6,2) - t3 * mrSges(7,1) - pkin(5) * t32 + t1 * mrSges(7,3) + qJD(6) * t71 + qJ(6) * t33 + t5 * mrSges(6,1) + t261; (m(7) * t133 - mrSges(7,2)) * t211 + (m(7) * t172 + t178 + t179) * t204 + (t133 * t154 + t158) * qJD(5) + t217; -mrSges(7,2) * t211 + qJD(5) * t158 - t153 * t246 + t217; t153; 0.2e1 * (m(7) * qJ(6) + mrSges(7,3)) * qJD(6); m(7) * t3 + t32; -m(7) * t273 - t201; (-m(7) * t246 - mrSges(7,2)) * t213; m(7) * t213; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
