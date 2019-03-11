% Calculate time derivative of joint inertia matrix for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:08:15
% EndTime: 2019-03-09 07:08:24
% DurationCPUTime: 4.08s
% Computational Cost: add. (12021->414), mult. (25773->609), div. (0->0), fcn. (27414->10), ass. (0->184)
t178 = sin(pkin(11));
t183 = sin(qJ(3));
t179 = cos(pkin(11));
t187 = cos(qJ(3));
t227 = t179 * t187;
t150 = -t183 * t178 + t227;
t151 = t178 * t187 + t183 * t179;
t182 = sin(qJ(4));
t186 = cos(qJ(4));
t117 = t150 * t182 + t151 * t186;
t181 = sin(qJ(5));
t224 = qJD(5) * t181;
t211 = t117 * t224;
t185 = cos(qJ(5));
t138 = t150 * qJD(3);
t139 = t151 * qJD(3);
t197 = t186 * t150 - t151 * t182;
t89 = qJD(4) * t197 + t138 * t186 - t139 * t182;
t234 = t185 * t89;
t192 = t211 - t234;
t223 = qJD(5) * t185;
t193 = t117 * t223 + t181 * t89;
t33 = mrSges(6,1) * t193 - mrSges(6,2) * t192;
t250 = pkin(7) + qJ(2);
t160 = t250 * t178;
t161 = t250 * t179;
t127 = -t183 * t160 + t187 * t161;
t109 = -t151 * qJD(2) - qJD(3) * t127;
t190 = -t138 * pkin(8) + t109;
t148 = t187 * t160;
t126 = -t161 * t183 - t148;
t111 = -pkin(8) * t151 + t126;
t112 = pkin(8) * t150 + t127;
t73 = t111 * t182 + t112 * t186;
t108 = -qJD(3) * t148 + qJD(2) * t227 + (-qJD(2) * t178 - qJD(3) * t161) * t183;
t97 = -pkin(8) * t139 + t108;
t41 = qJD(4) * t73 + t182 * t97 - t186 * t190;
t277 = m(6) * t41 + t33;
t269 = (t181 ^ 2 + t185 ^ 2) * t186;
t180 = sin(qJ(6));
t184 = cos(qJ(6));
t196 = t180 * t181 - t184 * t185;
t267 = qJD(5) + qJD(6);
t121 = t267 * t196;
t153 = t180 * t185 + t181 * t184;
t122 = t267 * t153;
t225 = -Ifges(7,5) * t121 - Ifges(7,6) * t122;
t276 = Ifges(6,5) * t223 + t225;
t270 = t186 * t111 - t112 * t182;
t40 = qJD(4) * t270 + t182 * t190 + t186 * t97;
t254 = pkin(3) * t139;
t90 = qJD(4) * t117 + t138 * t182 + t186 * t139;
t51 = pkin(4) * t90 - pkin(9) * t89 + t254;
t213 = -pkin(2) * t179 - pkin(1);
t132 = -pkin(3) * t150 + t213;
t74 = -pkin(4) * t197 - pkin(9) * t117 + t132;
t13 = t181 * t51 + t185 * t40 + t74 * t223 - t224 * t73;
t206 = -t181 * t40 + t185 * t51;
t67 = t185 * t73;
t48 = t181 * t74 + t67;
t14 = -qJD(5) * t48 + t206;
t275 = t13 * t185 - t14 * t181;
t273 = Ifges(6,5) * t234 + Ifges(6,3) * t90;
t169 = pkin(3) * t182 + pkin(9);
t249 = -pkin(10) - t169;
t145 = t249 * t181;
t173 = t185 * pkin(10);
t146 = t169 * t185 + t173;
t113 = t145 * t184 - t146 * t180;
t205 = qJD(5) * t249;
t245 = pkin(3) * qJD(4);
t216 = t186 * t245;
t130 = t181 * t205 + t185 * t216;
t131 = -t181 * t216 + t185 * t205;
t70 = qJD(6) * t113 + t130 * t184 + t131 * t180;
t114 = t145 * t180 + t146 * t184;
t71 = -qJD(6) * t114 - t130 * t180 + t131 * t184;
t272 = t71 * mrSges(7,1) - t70 * mrSges(7,2);
t258 = -pkin(10) - pkin(9);
t166 = t258 * t181;
t167 = pkin(9) * t185 + t173;
t128 = t166 * t184 - t167 * t180;
t212 = qJD(5) * t258;
t157 = t181 * t212;
t158 = t185 * t212;
t101 = qJD(6) * t128 + t157 * t184 + t158 * t180;
t129 = t166 * t180 + t167 * t184;
t102 = -qJD(6) * t129 - t157 * t180 + t158 * t184;
t271 = t102 * mrSges(7,1) - t101 * mrSges(7,2);
t79 = t196 * t117;
t47 = -t181 * t73 + t185 * t74;
t268 = -t181 * t47 + t185 * t48;
t202 = mrSges(6,1) * t181 + mrSges(6,2) * t185;
t154 = t202 * qJD(5);
t266 = 0.2e1 * m(6);
t265 = 2 * m(7);
t264 = -2 * mrSges(5,3);
t263 = 0.2e1 * t41;
t262 = -0.2e1 * t270;
t94 = t122 * mrSges(7,1) - t121 * mrSges(7,2);
t261 = 0.2e1 * t94;
t123 = mrSges(7,1) * t196 + mrSges(7,2) * t153;
t260 = 0.2e1 * t123;
t253 = pkin(3) * t186;
t252 = t41 * t270;
t248 = Ifges(6,4) * t181;
t247 = Ifges(6,4) * t185;
t246 = Ifges(6,6) * t181;
t244 = pkin(5) * qJD(6);
t242 = t122 * mrSges(7,3);
t237 = t182 * mrSges(5,1);
t233 = t186 * mrSges(5,2);
t231 = t117 * t181;
t230 = t117 * t185;
t163 = -mrSges(6,1) * t185 + mrSges(6,2) * t181;
t226 = t182 * t163;
t222 = qJD(6) * t180;
t221 = qJD(6) * t184;
t220 = 0.2e1 * t254;
t219 = 0.2e1 * mrSges(7,3);
t26 = -t117 * t122 - t196 * t89;
t27 = -t153 * t89 + t267 * t79;
t218 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t90;
t217 = mrSges(7,3) * t244;
t215 = pkin(5) * t224;
t214 = t184 * t121 * mrSges(7,3);
t171 = -pkin(5) * t185 - pkin(4);
t209 = t90 * mrSges(5,1) + t89 * mrSges(5,2);
t208 = -t224 / 0.2e1;
t207 = -(2 * Ifges(5,4)) - t246;
t203 = mrSges(6,3) * t269;
t201 = Ifges(6,1) * t185 - t248;
t200 = -Ifges(6,2) * t181 + t247;
t199 = Ifges(6,5) * t181 + Ifges(6,6) * t185;
t32 = -pkin(5) * t197 - pkin(10) * t230 + t47;
t39 = -pkin(10) * t231 + t48;
t16 = -t180 * t39 + t184 * t32;
t17 = t180 * t32 + t184 * t39;
t91 = mrSges(6,2) * t197 - mrSges(6,3) * t231;
t92 = -mrSges(6,1) * t197 - mrSges(6,3) * t230;
t198 = -t181 * t92 + t185 * t91;
t10 = -pkin(10) * t193 + t13;
t5 = -pkin(10) * t234 + pkin(5) * t90 + (-t67 + (pkin(10) * t117 - t74) * t181) * qJD(5) + t206;
t3 = qJD(6) * t16 + t10 * t184 + t180 * t5;
t4 = -qJD(6) * t17 - t10 * t180 + t184 * t5;
t195 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t218;
t194 = -t184 * t196 * t217 + (-pkin(5) * t242 + t153 * t217) * t180 + t276;
t124 = Ifges(7,4) * t153 - Ifges(7,2) * t196;
t125 = Ifges(7,1) * t153 - Ifges(7,4) * t196;
t155 = t200 * qJD(5);
t156 = t201 * qJD(5);
t164 = Ifges(6,2) * t185 + t248;
t165 = Ifges(6,1) * t181 + t247;
t95 = -Ifges(7,4) * t121 - Ifges(7,2) * t122;
t96 = -Ifges(7,1) * t121 - Ifges(7,4) * t122;
t191 = -t121 * t125 - t122 * t124 + t153 * t96 + t185 * t155 + t181 * t156 - t164 * t224 + t165 * t223 - t196 * t95;
t44 = mrSges(6,1) * t90 + mrSges(6,3) * t192;
t45 = -mrSges(6,2) * t90 - mrSges(6,3) * t193;
t189 = m(6) * (-t223 * t47 - t224 * t48 + t275) + t185 * t45 - t181 * t44 - t92 * t223 - t91 * t224;
t23 = pkin(5) * t193 + t41;
t30 = -Ifges(6,4) * t192 - Ifges(6,2) * t193 + Ifges(6,6) * t90;
t31 = -Ifges(6,1) * t192 - Ifges(6,4) * t193 + Ifges(6,5) * t90;
t78 = t153 * t117;
t42 = -Ifges(7,4) * t79 - Ifges(7,2) * t78 - Ifges(7,6) * t197;
t43 = -Ifges(7,1) * t79 - Ifges(7,4) * t78 - Ifges(7,5) * t197;
t55 = pkin(5) * t231 - t270;
t60 = -Ifges(6,6) * t197 + t117 * t200;
t61 = -Ifges(6,5) * t197 + t117 * t201;
t8 = Ifges(7,4) * t26 + Ifges(7,2) * t27 + Ifges(7,6) * t90;
t9 = Ifges(7,1) * t26 + Ifges(7,4) * t27 + Ifges(7,5) * t90;
t188 = t156 * t230 / 0.2e1 - t155 * t231 / 0.2e1 + t61 * t223 / 0.2e1 - t270 * t154 + (t163 - mrSges(5,1)) * t41 + (t117 * t208 + t234 / 0.2e1) * t165 - t17 * t242 + t60 * t208 + (Ifges(7,5) * t153 - Ifges(7,6) * t196 + t199) * t90 / 0.2e1 - t196 * t8 / 0.2e1 + (t121 * t16 - t153 * t4 - t196 * t3) * mrSges(7,3) + t185 * t30 / 0.2e1 + t181 * t31 / 0.2e1 + t153 * t9 / 0.2e1 + ((-t181 * t48 - t185 * t47) * qJD(5) + t275) * mrSges(6,3) - t193 * t164 / 0.2e1 - (-Ifges(6,6) * t224 + t276) * t197 / 0.2e1 - t40 * mrSges(5,2) + Ifges(5,5) * t89 - Ifges(5,6) * t90 + t55 * t94 - t78 * t95 / 0.2e1 - t79 * t96 / 0.2e1 - t121 * t43 / 0.2e1 - t122 * t42 / 0.2e1 + t23 * t123 + t27 * t124 / 0.2e1 + t26 * t125 / 0.2e1;
t170 = -pkin(4) - t253;
t162 = t171 - t253;
t159 = t182 * t245 + t215;
t144 = (-mrSges(7,1) * t180 - mrSges(7,2) * t184) * t244;
t136 = t138 * mrSges(4,2);
t83 = t202 * t117;
t59 = -mrSges(7,1) * t197 + mrSges(7,3) * t79;
t58 = mrSges(7,2) * t197 - mrSges(7,3) * t78;
t52 = mrSges(7,1) * t78 - mrSges(7,2) * t79;
t19 = -mrSges(7,2) * t90 + mrSges(7,3) * t27;
t18 = mrSges(7,1) * t90 - mrSges(7,3) * t26;
t11 = -mrSges(7,1) * t27 + mrSges(7,2) * t26;
t1 = [0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t178 ^ 2 + t179 ^ 2) * qJD(2) + (mrSges(5,3) * t262 - t181 * t60 + t185 * t61) * t89 + 0.2e1 * t132 * t209 + t33 * t262 + t83 * t263 + (t16 * t4 + t17 * t3 + t23 * t55) * t265 + 0.2e1 * t213 * (t139 * mrSges(4,1) + t136) + 0.2e1 * (t138 * t150 - t139 * t151) * Ifges(4,4) + 0.2e1 * (t108 * t150 - t109 * t151 - t126 * t138 - t127 * t139) * mrSges(4,3) + 0.2e1 * m(5) * (t132 * t254 + t40 * t73 - t252) + (t13 * t48 + t14 * t47 - t252) * t266 - (mrSges(5,1) * t220 + t207 * t89 + t40 * t264 + t218 + t273) * t197 + (mrSges(5,2) * t220 + mrSges(5,3) * t263 + 0.2e1 * Ifges(5,1) * t89 - t181 * t30 + t185 * t31 + (-t181 * t61 - t185 * t60 + t197 * t199) * qJD(5)) * t117 + ((Ifges(6,5) * t185 + t207) * t117 - ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t197 + t73 * t264 - Ifges(7,5) * t79 - Ifges(7,6) * t78) * t90 - 0.2e1 * t150 * Ifges(4,2) * t139 + 0.2e1 * t151 * t138 * Ifges(4,1) + 0.2e1 * m(4) * (t108 * t127 + t109 * t126) + 0.2e1 * t16 * t18 + 0.2e1 * t17 * t19 + t27 * t42 + t26 * t43 + 0.2e1 * t47 * t44 + 0.2e1 * t48 * t45 + 0.2e1 * t23 * t52 + 0.2e1 * t55 * t11 + 0.2e1 * t3 * t58 + 0.2e1 * t4 * t59 - t78 * t8 - t79 * t9 + 0.2e1 * t13 * t91 + 0.2e1 * t14 * t92; -t121 * t58 - t122 * t59 - t196 * t18 + t153 * t19 + t181 * t45 + t185 * t44 + t136 - (-m(5) * pkin(3) - mrSges(4,1)) * t139 + t198 * qJD(5) + m(7) * (-t121 * t17 - t122 * t16 + t153 * t3 - t196 * t4) + m(6) * (qJD(5) * t268 + t13 * t181 + t14 * t185) + t209; (-t121 * t153 + t122 * t196) * t265; m(7) * (t113 * t4 + t114 * t3 + t159 * t55 + t16 * t71 + t162 * t23 + t17 * t70) + (m(5) * (t182 * t40 - t186 * t41) + (-t182 * t90 - t186 * t89) * mrSges(5,3) + ((m(5) * t73 + m(6) * t268 + mrSges(5,3) * t197 + t198) * t186 + (t117 * mrSges(5,3) + t83 - (m(6) + m(5)) * t270) * t182) * qJD(4)) * pkin(3) + t189 * t169 + t159 * t52 + t162 * t11 + t188 + t70 * t58 + t71 * t59 - t108 * mrSges(4,2) + t109 * mrSges(4,1) + t113 * t18 + t114 * t19 + Ifges(4,5) * t138 - Ifges(4,6) * t139 + t277 * t170; m(7) * (-t113 * t122 - t114 * t121 + t153 * t70 - t196 * t71); t159 * t260 + t162 * t261 + (t113 * t71 + t114 * t70 + t159 * t162) * t265 + 0.2e1 * t170 * t154 + (t113 * t121 - t114 * t122 - t71 * t153 - t196 * t70) * t219 + (-0.2e1 * t233 - 0.2e1 * t237 + 0.2e1 * t226 + (t169 * t269 + t170 * t182) * t266 + 0.2e1 * t203) * t245 + t191; m(7) * (t101 * t17 + t102 * t16 + t128 * t4 + t129 * t3 + t171 * t23 + t215 * t55) + t189 * pkin(9) + t171 * t11 + t52 * t215 + t188 + t101 * t58 + t102 * t59 + t128 * t18 + t129 * t19 - t277 * pkin(4); m(7) * (t101 * t153 - t102 * t196 - t121 * t129 - t122 * t128); m(7) * (t101 * t114 + t102 * t113 + t128 * t71 + t129 * t70 + t159 * t171 + t162 * t215) + (t171 + t162) * t94 + (t170 - pkin(4)) * t154 + (t159 + t215) * t123 + (m(6) * (-pkin(4) * t182 + pkin(9) * t269) - t233 + t226 - t237 + t203) * t245 + ((-t102 - t71) * t153 - (t101 + t70) * t196 - (t114 + t129) * t122 - (-t113 - t128) * t121) * mrSges(7,3) + t191; t215 * t260 + t171 * t261 - 0.2e1 * pkin(4) * t154 + (t101 * t129 + t102 * t128 + t171 * t215) * t265 + (-t101 * t196 - t102 * t153 + t128 * t121 - t129 * t122) * t219 + t191; -Ifges(6,5) * t211 + t14 * mrSges(6,1) - t13 * mrSges(6,2) - t193 * Ifges(6,6) + (m(7) * (-t16 * t222 + t17 * t221 + t180 * t3 + t184 * t4) + t58 * t221 + t180 * t19 - t59 * t222 + t184 * t18) * pkin(5) + t195 + t273; -t154 + m(7) * (-t121 * t180 - t122 * t184 + (t153 * t184 + t180 * t196) * qJD(6)) * pkin(5) - t94; -t202 * t216 + (t163 * t169 - t246) * qJD(5) + (m(7) * (-t113 * t222 + t114 * t221 + t180 * t70 + t184 * t71) + t214) * pkin(5) + t194 + t272; (pkin(9) * t163 - t246) * qJD(5) + (m(7) * (t101 * t180 + t102 * t184 - t128 * t222 + t129 * t221) + t214) * pkin(5) + t194 + t271; 0.2e1 * t144; t195; -t94; t225 + t272; t225 + t271; t144; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
