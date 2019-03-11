% Calculate time derivative of joint inertia matrix for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:16:59
% EndTime: 2019-03-09 07:17:07
% DurationCPUTime: 3.69s
% Computational Cost: add. (7458->333), mult. (14643->502), div. (0->0), fcn. (14097->8), ass. (0->168)
t233 = qJD(3) + qJD(4);
t123 = sin(qJ(5));
t127 = cos(qJ(5));
t124 = sin(qJ(4));
t125 = sin(qJ(3));
t128 = cos(qJ(4));
t129 = cos(qJ(3));
t95 = -t124 * t129 - t128 * t125;
t180 = t128 * t129;
t183 = t124 * t125;
t96 = t180 - t183;
t149 = t123 * t96 - t127 * t95;
t178 = qJD(4) * t124;
t74 = -qJD(3) * t183 - t125 * t178 + t180 * t233;
t75 = t233 * t95;
t36 = t149 * qJD(5) + t123 * t74 - t127 * t75;
t130 = -pkin(1) - pkin(7);
t211 = pkin(8) - t130;
t168 = t129 * t211;
t160 = t128 * t168;
t238 = t211 * t125;
t76 = t124 * t238 - t160;
t139 = -t96 * pkin(9) + t76;
t244 = t124 * t168 + t128 * t238;
t57 = pkin(9) * t95 - t244;
t40 = t123 * t57 - t127 * t139;
t250 = t36 * t40;
t69 = t123 * t95 + t127 * t96;
t219 = t36 * t69;
t122 = sin(qJ(6));
t120 = t122 ^ 2;
t126 = cos(qJ(6));
t121 = t126 ^ 2;
t179 = t120 + t121;
t136 = qJD(5) * t69 + t123 * t75 + t127 * t74;
t249 = -t123 * t136 + t127 * t36;
t116 = pkin(3) * t128 + pkin(4);
t176 = qJD(5) * t127;
t177 = qJD(5) * t123;
t182 = t124 * t127;
t63 = t116 * t177 + (t124 * t176 + (t123 * t128 + t182) * qJD(4)) * pkin(3);
t214 = t63 * t69;
t184 = t123 * t124;
t62 = t116 * t176 + (-t124 * t177 + (t127 * t128 - t184) * qJD(4)) * pkin(3);
t215 = t62 * t149;
t84 = -pkin(3) * t184 + t116 * t127;
t85 = pkin(3) * t182 + t123 * t116;
t248 = -t136 * t85 + t36 * t84 + t214 - t215;
t208 = mrSges(7,1) * t126;
t101 = mrSges(7,2) * t122 - t208;
t247 = t101 - mrSges(6,1);
t174 = qJD(6) * t126;
t199 = t122 * t36;
t146 = t69 * t174 - t199;
t245 = qJD(3) * t129;
t41 = t123 * t139 + t127 * t57;
t243 = t136 * t41;
t47 = qJD(3) * t76 - qJD(4) * t160 + t238 * t178;
t241 = -pkin(9) * t74 + qJD(5) * t139 + t47;
t239 = t136 * t149;
t237 = t179 * t127;
t114 = pkin(4) * t123 + pkin(10);
t236 = t179 * t114;
t111 = t125 * pkin(3) + qJ(2);
t78 = -pkin(4) * t95 + t111;
t42 = pkin(5) * t149 - pkin(10) * t69 + t78;
t16 = -t122 * t41 + t126 * t42;
t17 = t122 * t42 + t126 * t41;
t235 = -t122 * t16 + t126 * t17;
t234 = -mrSges(4,1) * t125 - mrSges(4,2) * t129;
t48 = t233 * t244;
t232 = t48 * mrSges(5,1) - t47 * mrSges(5,2) + Ifges(5,5) * t75 - Ifges(5,6) * t74;
t172 = pkin(4) * t176;
t161 = mrSges(7,3) * t172;
t104 = t120 * t161;
t105 = t121 * t161;
t115 = -pkin(4) * t127 - pkin(5);
t155 = mrSges(7,1) * t122 + mrSges(7,2) * t126;
t97 = t155 * qJD(6);
t79 = t115 * t97;
t91 = pkin(4) * t101 * t177;
t231 = t104 + t105 + t79 + t91;
t230 = -t124 * t74 - t128 * t75 + (t124 * t96 + t128 * t95) * qJD(4);
t229 = 2 * m(6);
t228 = 2 * m(7);
t133 = -t75 * pkin(9) + t48;
t11 = t123 * t241 - t127 * t133 + t57 * t176;
t227 = 0.2e1 * t11;
t107 = pkin(3) * t245 + qJD(2);
t60 = pkin(4) * t74 + t107;
t226 = 0.2e1 * t60;
t225 = m(6) / 0.2e1;
t224 = pkin(5) * t97;
t222 = t11 * t40;
t10 = t123 * t133 + t241 * t127 - t57 * t177;
t15 = pkin(5) * t136 + pkin(10) * t36 + t60;
t191 = qJD(6) * t16;
t2 = t10 * t126 + t122 * t15 + t191;
t221 = t126 * t2;
t190 = qJD(6) * t17;
t3 = -t10 * t122 + t126 * t15 - t190;
t220 = t3 * t122;
t218 = t40 * t63;
t216 = t62 * mrSges(6,2);
t213 = t74 * t95;
t212 = t75 * t96;
t193 = t126 * t36;
t210 = -Ifges(7,5) * t193 + Ifges(7,3) * t136;
t206 = mrSges(7,3) * t120;
t205 = mrSges(7,3) * t121;
t204 = Ifges(7,4) * t122;
t203 = Ifges(7,4) * t126;
t202 = Ifges(7,6) * t122;
t201 = pkin(4) * qJD(5);
t198 = t122 * t69;
t196 = t69 * t123;
t192 = t126 * t69;
t81 = pkin(10) + t85;
t189 = qJD(6) * t81;
t175 = qJD(6) * t122;
t173 = 2 * mrSges(5,3);
t171 = t69 * t175;
t145 = t171 + t193;
t12 = t146 * mrSges(7,1) - t145 * mrSges(7,2);
t169 = m(7) * t11 + t12;
t166 = -t175 / 0.2e1;
t165 = -(2 * Ifges(6,4)) - t202;
t163 = t62 * t179;
t162 = t179 * t136;
t159 = -t11 * t69 + t250;
t157 = -t212 + t213;
t156 = -t123 * mrSges(6,1) - t127 * mrSges(6,2);
t154 = Ifges(7,1) * t126 - t204;
t153 = -Ifges(7,2) * t122 + t203;
t152 = Ifges(7,5) * t122 + Ifges(7,6) * t126;
t13 = mrSges(7,1) * t136 + t145 * mrSges(7,3);
t14 = -mrSges(7,2) * t136 - t146 * mrSges(7,3);
t151 = -t122 * t13 + t126 * t14;
t44 = -mrSges(7,2) * t149 - mrSges(7,3) * t198;
t45 = mrSges(7,1) * t149 - mrSges(7,3) * t192;
t150 = -t122 * t45 + t126 * t44;
t102 = Ifges(7,2) * t126 + t204;
t103 = Ifges(7,1) * t122 + t203;
t98 = t153 * qJD(6);
t99 = t154 * qJD(6);
t144 = -t102 * t175 + t103 * t174 + t122 * t99 + t126 * t98;
t143 = -t69 * t97 + t247 * t36 + (-mrSges(6,2) + t205 + t206) * t136;
t142 = (-mrSges(5,1) * t124 - mrSges(5,2) * t128) * qJD(4) * pkin(3);
t141 = t244 * t74 + t47 * t95 - t48 * t96 - t75 * t76;
t140 = -t220 + (-t122 * t17 - t126 * t16) * qJD(6);
t138 = t75 * mrSges(5,1) - t74 * mrSges(5,2) + t143;
t55 = t63 * t101;
t58 = t62 * t206;
t59 = t62 * t205;
t61 = t63 * mrSges(6,1);
t80 = -pkin(5) - t84;
t73 = t80 * t97;
t135 = t144 + t55 + t58 + t59 - t61 + t73 - t216;
t117 = Ifges(7,5) * t174;
t24 = Ifges(7,6) * t149 + t153 * t69;
t25 = Ifges(7,5) * t149 + t154 * t69;
t6 = -t145 * Ifges(7,4) - t146 * Ifges(7,2) + Ifges(7,6) * t136;
t7 = -t145 * Ifges(7,1) - t146 * Ifges(7,4) + Ifges(7,5) * t136;
t134 = -t10 * mrSges(6,2) + mrSges(7,3) * t221 + t24 * t166 + t25 * t174 / 0.2e1 + t40 * t97 - Ifges(6,5) * t36 + t122 * t7 / 0.2e1 - t98 * t198 / 0.2e1 + t126 * t6 / 0.2e1 + t99 * t192 / 0.2e1 + t149 * (-Ifges(7,6) * t175 + t117) / 0.2e1 - t146 * t102 / 0.2e1 + (t152 / 0.2e1 - Ifges(6,6)) * t136 + t247 * t11 + (-t193 / 0.2e1 + t69 * t166) * t103;
t132 = -t45 * t174 - t44 * t175 + m(7) * (-t16 * t174 - t17 * t175 - t220 + t221) + t151;
t131 = t140 * mrSges(7,3) + t134;
t43 = t155 * t69;
t1 = [-t25 * t193 + t24 * t199 + 0.2e1 * t16 * t13 + 0.2e1 * t17 * t14 + 0.2e1 * t40 * t12 + t43 * t227 + 0.2e1 * t2 * t44 + 0.2e1 * t3 * t45 + 0.2e1 * t78 * (mrSges(6,1) * t136 - mrSges(6,2) * t36) - 0.2e1 * Ifges(5,2) * t213 + 0.2e1 * Ifges(5,1) * t212 + 0.2e1 * t107 * (-mrSges(5,1) * t95 + mrSges(5,2) * t96) + 0.2e1 * t111 * (mrSges(5,1) * t74 + mrSges(5,2) * t75) + 0.2e1 * (-t243 - t250) * mrSges(6,3) + 0.2e1 * m(5) * (t107 * t111 - t244 * t47 + t48 * t76) + (t10 * t41 + t60 * t78 + t222) * t229 + (t16 * t3 + t17 * t2 + t222) * t228 + 0.2e1 * (-t74 * t96 + t75 * t95) * Ifges(5,4) + t141 * t173 + (mrSges(6,1) * t226 - 0.2e1 * mrSges(6,3) * t10 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t136 - t165 * t36 + t210) * t149 + (mrSges(6,2) * t226 + mrSges(6,3) * t227 - 0.2e1 * Ifges(6,1) * t36 - t122 * t6 + t126 * t7 + (Ifges(7,5) * t126 + t165) * t136 + (-t122 * t25 - t126 * t24 - t149 * t152) * qJD(6)) * t69 + 0.2e1 * (qJ(2) * (mrSges(4,1) * t129 - mrSges(4,2) * t125) + (t125 ^ 2 - t129 ^ 2) * Ifges(4,4)) * qJD(3) + 0.2e1 * (mrSges(3,3) + (m(3) + m(4)) * qJ(2) - t234) * qJD(2) + 0.2e1 * (-Ifges(4,1) + Ifges(4,2)) * t125 * t245; -t69 * t12 + t36 * t43 + t150 * t136 + t157 * t173 + ((-t122 * t44 - t126 * t45) * qJD(6) + t151) * t149 + m(7) * (t235 * t136 + (t140 + t221) * t149 + t159) + m(6) * (t10 * t149 + t159 + t243) - m(5) * t141 + (0.2e1 * t219 - 0.2e1 * t239) * mrSges(6,3); 0.2e1 * m(7) * (t149 * t162 - t219) + 0.2e1 * m(6) * (-t219 + t239) - 0.2e1 * m(5) * t157; t134 + (-t44 * t189 + m(7) * (-t16 * t62 - t17 * t189 - t3 * t81) - t81 * t13 - t62 * t45 + (-t3 - t190) * mrSges(7,3)) * t122 + (-mrSges(7,3) * t191 - t45 * t189 + m(7) * (-t16 * t189 + t17 * t62 + t2 * t81) + t81 * t14 + t62 * t44) * t126 + (m(5) * (t124 * t47 + t128 * t48 + (-t124 * t76 - t128 * t244) * qJD(4)) + t230 * mrSges(5,3)) * pkin(3) + ((-mrSges(4,2) * t130 - Ifges(4,6)) * t129 + (-mrSges(4,1) * t130 - Ifges(4,5)) * t125) * qJD(3) + m(7) * (t11 * t80 + t218) + m(6) * (t10 * t85 - t11 * t84 + t41 * t62 + t218) + t248 * mrSges(6,3) + t63 * t43 + t80 * t12 + t232; t234 * qJD(3) + m(7) * (t36 * t80 - t214 + t179 * (t136 * t81 + t215)) - m(6) * t248 - m(5) * t230 * pkin(3) + t138; -0.2e1 * t216 + 0.2e1 * t55 + 0.2e1 * t58 + 0.2e1 * t59 - 0.2e1 * t61 + 0.2e1 * t73 + 0.2e1 * t142 + (t81 * t163 + t63 * t80) * t228 + (t62 * t85 - t63 * t84) * t229 + t144; t131 + t132 * t114 + t169 * t115 + (m(6) * (t10 * t123 - t11 * t127) + t249 * mrSges(6,3) + ((m(6) * t41 + m(7) * t235 - t149 * mrSges(6,3) + t150) * t127 + (t69 * mrSges(6,3) + t43 + (m(7) + m(6)) * t40) * t123) * qJD(5)) * pkin(4) + t232; m(7) * (t115 * t36 + t136 * t236) + 0.2e1 * (-t249 * t225 + (m(7) * (t149 * t237 - t196) / 0.2e1 + (t127 * t149 - t196) * t225) * qJD(5)) * pkin(4) + t138; m(7) * (t115 * t63 + t236 * t62) + (m(6) * (t123 * t62 - t127 * t63) + (m(7) * (t123 * t80 + t237 * t81) + m(6) * (-t123 * t84 + t127 * t85) + t156) * qJD(5)) * pkin(4) + t135 + t142 + t231; 0.2e1 * t104 + 0.2e1 * t105 + 0.2e1 * t79 + 0.2e1 * t91 + 0.2e1 * (m(7) * (t114 * t237 + t115 * t123) + t156) * t201 + t144; -t169 * pkin(5) + t132 * pkin(10) + t131; m(7) * (-pkin(5) * t36 + pkin(10) * t162) + t143; m(7) * (-pkin(5) * t63 + pkin(10) * t163) - t224 + t135; -t224 + (m(7) * (-pkin(5) * t123 + pkin(10) * t237) + t156) * t201 + t144 + t231; t144 - 0.2e1 * t224; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t171 - t146 * Ifges(7,6) + t210; (-t126 * t136 + t149 * t175) * mrSges(7,2) + (-t122 * t136 - t149 * t174) * mrSges(7,1); t117 - t155 * t62 + (-t81 * t208 + (mrSges(7,2) * t81 - Ifges(7,6)) * t122) * qJD(6); t117 - t155 * t172 + (t101 * t114 - t202) * qJD(6); t117 + (pkin(10) * t101 - t202) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
