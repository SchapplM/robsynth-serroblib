% Calculate time derivative of joint inertia matrix for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:49:35
% EndTime: 2019-03-09 10:49:43
% DurationCPUTime: 4.13s
% Computational Cost: add. (5039->470), mult. (11887->682), div. (0->0), fcn. (10676->8), ass. (0->179)
t238 = Ifges(6,4) + Ifges(5,5);
t237 = -Ifges(6,2) - Ifges(5,3);
t236 = Ifges(6,6) - Ifges(5,6);
t178 = sin(pkin(10));
t179 = cos(pkin(10));
t235 = mrSges(4,1) * t178 + mrSges(4,2) * t179;
t182 = sin(qJ(2));
t205 = qJD(2) * t182;
t223 = pkin(8) + qJ(3);
t161 = t223 * t178;
t162 = t223 * t179;
t181 = sin(qJ(4));
t184 = cos(qJ(4));
t117 = -t161 * t181 + t162 * t184;
t212 = t178 * t181;
t152 = -t179 * t184 + t212;
t141 = t152 * qJD(4);
t153 = t178 * t184 + t179 * t181;
t142 = t153 * qJD(4);
t234 = -t141 * t238 + t142 * t236;
t180 = sin(qJ(6));
t183 = cos(qJ(6));
t233 = -mrSges(7,1) * t180 - mrSges(7,2) * t183;
t200 = qJD(4) * t182;
t185 = cos(qJ(2));
t204 = qJD(2) * t185;
t86 = -t152 * t204 - t153 * t200;
t199 = qJD(4) * t184;
t209 = t179 * t182;
t87 = t153 * t204 + t199 * t209 - t200 * t212;
t232 = t205 * t237 - t236 * t87 - t238 * t86;
t231 = 2 * m(4);
t230 = 2 * m(5);
t229 = 2 * m(6);
t228 = 2 * m(7);
t227 = -2 * pkin(1);
t226 = 2 * pkin(7);
t186 = -pkin(4) - pkin(5);
t225 = t179 / 0.2e1;
t131 = t153 * t182;
t132 = t152 * t182;
t72 = t131 * t183 + t132 * t180;
t24 = qJD(6) * t72 + t180 * t87 + t183 * t86;
t73 = t131 * t180 - t132 * t183;
t25 = -qJD(6) * t73 - t180 * t86 + t183 * t87;
t222 = Ifges(7,5) * t24 + Ifges(7,6) * t25;
t102 = t152 * t183 - t153 * t180;
t42 = qJD(6) * t102 - t141 * t183 + t142 * t180;
t103 = t152 * t180 + t153 * t183;
t43 = -qJD(6) * t103 + t141 * t180 + t142 * t183;
t221 = Ifges(7,5) * t42 + Ifges(7,6) * t43;
t216 = Ifges(4,4) * t178;
t215 = Ifges(4,4) * t179;
t157 = -qJ(5) * t180 + t183 * t186;
t129 = qJD(5) * t183 + qJD(6) * t157;
t214 = t129 * mrSges(7,2);
t158 = qJ(5) * t183 + t180 * t186;
t130 = -qJD(5) * t180 - qJD(6) * t158;
t213 = t130 * mrSges(7,1);
t159 = -pkin(2) * t185 - qJ(3) * t182 - pkin(1);
t148 = t179 * t159;
t105 = -pkin(8) * t209 + t148 + (-pkin(7) * t178 - pkin(3)) * t185;
t208 = t179 * t185;
t125 = pkin(7) * t208 + t159 * t178;
t211 = t178 * t182;
t115 = -pkin(8) * t211 + t125;
t52 = t105 * t181 + t115 * t184;
t210 = t178 * t185;
t133 = t235 * t204;
t140 = -qJD(3) * t182 + (pkin(2) * t182 - qJ(3) * t185) * qJD(2);
t198 = pkin(7) * t205;
t113 = t140 * t179 + t178 * t198;
t156 = pkin(3) * t211 + pkin(7) * t182;
t203 = qJD(3) * t178;
t202 = qJD(3) * t179;
t201 = qJD(4) * t181;
t170 = -pkin(3) * t179 - pkin(2);
t39 = t87 * mrSges(5,1) + mrSges(5,2) * t86;
t38 = mrSges(6,1) * t87 - t86 * mrSges(6,3);
t116 = t161 * t184 + t162 * t181;
t66 = -t161 * t199 + t184 * t202 + (-qJD(4) * t162 - t203) * t181;
t67 = qJD(3) * t153 + qJD(4) * t117;
t197 = t116 * t67 + t117 * t66;
t96 = t142 * mrSges(5,1) - mrSges(5,2) * t141;
t95 = mrSges(6,1) * t142 + t141 * mrSges(6,3);
t51 = t105 * t184 - t115 * t181;
t63 = -mrSges(6,1) * t205 + mrSges(6,2) * t86;
t44 = -qJ(5) * t185 + t52;
t5 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t17 = -t43 * mrSges(7,1) + t42 * mrSges(7,2);
t195 = -qJ(5) * t132 - t156;
t46 = pkin(4) * t185 - t51;
t194 = t179 * Ifges(4,1) - t216;
t193 = -t178 * Ifges(4,2) + t215;
t192 = -Ifges(4,5) * t179 + Ifges(4,6) * t178;
t27 = pkin(5) * t185 + pkin(9) * t132 + t46;
t30 = pkin(9) * t131 + t44;
t10 = -t180 * t30 + t183 * t27;
t11 = t180 * t27 + t183 * t30;
t84 = -pkin(9) * t153 + t116;
t85 = pkin(9) * t152 + t117;
t36 = -t180 * t85 + t183 * t84;
t37 = t180 * t84 + t183 * t85;
t74 = (pkin(3) * t182 - pkin(8) * t208) * qJD(2) + t113;
t126 = t178 * t140;
t91 = t126 + (-pkin(7) * t209 - pkin(8) * t210) * qJD(2);
t16 = -t105 * t201 - t115 * t199 - t181 * t91 + t184 * t74;
t191 = -qJ(5) * t141 + qJD(5) * t153;
t190 = qJ(5) * t153 - t170;
t145 = (pkin(3) * t178 + pkin(7)) * t204;
t54 = pkin(9) * t142 + t66;
t55 = pkin(9) * t141 + t67;
t6 = qJD(6) * t36 + t180 * t55 + t183 * t54;
t7 = -qJD(6) * t37 - t180 * t54 + t183 * t55;
t189 = t7 * mrSges(7,1) - t6 * mrSges(7,2) + t221;
t15 = t105 * t199 - t115 * t201 + t181 * t74 + t184 * t91;
t26 = pkin(4) * t87 - qJ(5) * t86 + qJD(5) * t132 + t145;
t8 = -pkin(9) * t86 + t186 * t205 - t16;
t12 = qJ(5) * t205 - qJD(5) * t185 + t15;
t9 = pkin(9) * t87 + t12;
t1 = qJD(6) * t10 + t180 * t8 + t183 * t9;
t2 = -qJD(6) * t11 - t180 * t9 + t183 * t8;
t188 = -t2 * mrSges(7,1) + t1 * mrSges(7,2) + Ifges(7,3) * t205 - t222;
t155 = -mrSges(4,1) * t185 - mrSges(4,3) * t209;
t154 = mrSges(4,2) * t185 - mrSges(4,3) * t211;
t144 = (mrSges(4,1) * t182 - mrSges(4,3) * t208) * qJD(2);
t143 = (-mrSges(4,2) * t182 - mrSges(4,3) * t210) * qJD(2);
t124 = -pkin(7) * t210 + t148;
t123 = (t182 * Ifges(4,5) + t185 * t194) * qJD(2);
t122 = (t182 * Ifges(4,6) + t185 * t193) * qJD(2);
t121 = mrSges(6,1) * t185 - mrSges(6,2) * t132;
t120 = -mrSges(5,1) * t185 + mrSges(5,3) * t132;
t119 = mrSges(5,2) * t185 - mrSges(5,3) * t131;
t118 = -mrSges(6,2) * t131 - mrSges(6,3) * t185;
t114 = -t179 * t198 + t126;
t112 = Ifges(5,1) * t153 - Ifges(5,4) * t152;
t111 = Ifges(6,1) * t153 + Ifges(6,5) * t152;
t110 = Ifges(5,4) * t153 - Ifges(5,2) * t152;
t109 = Ifges(6,5) * t153 + Ifges(6,3) * t152;
t108 = mrSges(6,1) * t152 - mrSges(6,3) * t153;
t101 = pkin(4) * t152 - t190;
t100 = -Ifges(5,1) * t141 - Ifges(5,4) * t142;
t99 = -Ifges(6,1) * t141 + Ifges(6,5) * t142;
t98 = -Ifges(5,4) * t141 - Ifges(5,2) * t142;
t97 = -Ifges(6,5) * t141 + Ifges(6,3) * t142;
t90 = mrSges(6,1) * t131 + mrSges(6,3) * t132;
t71 = -Ifges(5,1) * t132 - Ifges(5,4) * t131 - Ifges(5,5) * t185;
t70 = -Ifges(6,1) * t132 - Ifges(6,4) * t185 + Ifges(6,5) * t131;
t69 = -Ifges(5,4) * t132 - Ifges(5,2) * t131 - Ifges(5,6) * t185;
t68 = -Ifges(6,5) * t132 - Ifges(6,6) * t185 + Ifges(6,3) * t131;
t64 = -mrSges(5,2) * t205 - mrSges(5,3) * t87;
t62 = mrSges(5,1) * t205 - mrSges(5,3) * t86;
t61 = -mrSges(6,2) * t87 + mrSges(6,3) * t205;
t60 = t152 * t186 + t190;
t59 = mrSges(7,1) * t185 - mrSges(7,3) * t73;
t58 = -mrSges(7,2) * t185 + mrSges(7,3) * t72;
t57 = pkin(4) * t142 - t191;
t56 = pkin(4) * t131 - t195;
t53 = t142 * t186 + t191;
t50 = Ifges(7,1) * t103 + Ifges(7,4) * t102;
t49 = Ifges(7,4) * t103 + Ifges(7,2) * t102;
t48 = -mrSges(7,1) * t102 + mrSges(7,2) * t103;
t47 = t131 * t186 + t195;
t35 = -mrSges(7,1) * t72 + mrSges(7,2) * t73;
t34 = Ifges(5,1) * t86 - Ifges(5,4) * t87 + Ifges(5,5) * t205;
t33 = Ifges(6,1) * t86 + Ifges(6,4) * t205 + Ifges(6,5) * t87;
t32 = Ifges(5,4) * t86 - Ifges(5,2) * t87 + Ifges(5,6) * t205;
t31 = Ifges(6,5) * t86 + Ifges(6,6) * t205 + Ifges(6,3) * t87;
t29 = Ifges(7,1) * t73 + Ifges(7,4) * t72 + Ifges(7,5) * t185;
t28 = Ifges(7,4) * t73 + Ifges(7,2) * t72 + Ifges(7,6) * t185;
t21 = mrSges(7,2) * t205 + mrSges(7,3) * t25;
t20 = -mrSges(7,1) * t205 - mrSges(7,3) * t24;
t19 = Ifges(7,1) * t42 + Ifges(7,4) * t43;
t18 = Ifges(7,4) * t42 + Ifges(7,2) * t43;
t14 = pkin(5) * t87 + t26;
t13 = -pkin(4) * t205 - t16;
t4 = Ifges(7,1) * t24 + Ifges(7,4) * t25 - Ifges(7,5) * t205;
t3 = Ifges(7,4) * t24 + Ifges(7,2) * t25 - Ifges(7,6) * t205;
t22 = [(t1 * t11 + t10 * t2 - t14 * t47) * t228 + (t12 * t44 + t13 * t46 + t26 * t56) * t229 + (t145 * t156 + t15 * t52 + t16 * t51) * t230 + (t113 * t124 + t114 * t125) * t231 + (t68 - t69) * t87 + (t70 + t71) * t86 + 0.2e1 * t156 * t39 + 0.2e1 * t114 * t154 + 0.2e1 * t113 * t155 + 0.2e1 * t125 * t143 + 0.2e1 * t124 * t144 - t131 * t32 - t132 * t33 - t132 * t34 + t131 * t31 + 0.2e1 * t12 * t118 + 0.2e1 * t15 * t119 + 0.2e1 * t16 * t120 + 0.2e1 * t13 * t121 + 0.2e1 * t26 * t90 + t72 * t3 + t73 * t4 + 0.2e1 * t2 * t59 + 0.2e1 * t44 * t61 + 0.2e1 * t51 * t62 + 0.2e1 * t46 * t63 + 0.2e1 * t52 * t64 + 0.2e1 * t56 * t38 + 0.2e1 * t1 * t58 + 0.2e1 * t47 * t5 - 0.2e1 * t14 * t35 + t25 * t28 + t24 * t29 + 0.2e1 * t10 * t20 + 0.2e1 * t11 * t21 + (-t122 * t178 + t123 * t179 + t133 * t226) * t182 + (t222 + t232) * t185 + 0.2e1 * t145 * (mrSges(5,1) * t131 - mrSges(5,2) * t132) + (((mrSges(3,2) * t227) + 0.2e1 * (Ifges(3,4) + t192) * t185) * t185 + ((mrSges(3,1) * t227) - Ifges(7,5) * t73 - Ifges(7,6) * t72 + (-0.2e1 * Ifges(3,4) - t192) * t182 - t238 * t132 + t236 * t131 + ((pkin(7) ^ 2 * t231) - t178 * t193 + t179 * t194 + t226 * t235 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - (2 * Ifges(4,3)) - 0.2e1 * Ifges(7,3) + t237) * t185) * t182) * qJD(2); (t122 / 0.2e1 + qJ(3) * t143 + qJD(3) * t154 + t114 * mrSges(4,3)) * t179 + (t123 / 0.2e1 - qJ(3) * t144 - qJD(3) * t155 - t113 * mrSges(4,3)) * t178 + (t1 * t102 - t10 * t42 - t103 * t2 + t11 * t43) * mrSges(7,3) + (t61 + t64) * t117 + (t63 - t62) * t116 + m(5) * (-t116 * t16 + t117 * t15 + t145 * t170 - t51 * t67 + t52 * t66) + m(6) * (t101 * t26 + t116 * t13 + t117 * t12 + t44 * t66 + t46 * t67 + t56 * t57) + m(7) * (t1 * t37 + t10 * t7 + t11 * t6 - t14 * t60 + t2 * t36 + t47 * t53) + (t109 / 0.2e1 - t110 / 0.2e1) * t87 + (-t120 + t121) * t67 + (t118 + t119) * t66 + (t111 / 0.2e1 + t112 / 0.2e1) * t86 + t156 * t96 + t170 * t39 - pkin(2) * t133 + t101 * t38 + t102 * t3 / 0.2e1 + t103 * t4 / 0.2e1 + t26 * t108 + t57 * t90 + t56 * t95 + t72 * t18 / 0.2e1 + t73 * t19 / 0.2e1 + t7 * t59 + t60 * t5 + t24 * t50 / 0.2e1 + t53 * t35 + t6 * t58 + t42 * t29 / 0.2e1 + t43 * t28 / 0.2e1 + t47 * t17 - t14 * t48 + t25 * t49 / 0.2e1 + t36 * t20 + t37 * t21 + ((pkin(7) * mrSges(3,2)) + Ifges(4,5) * t178 / 0.2e1 + Ifges(4,6) * t225 - Ifges(3,6) - Ifges(7,5) * t103 / 0.2e1 - Ifges(7,6) * t102 / 0.2e1 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t153 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t152) * t205 + (t141 * t51 - t142 * t52 - t15 * t152 - t153 * t16) * mrSges(5,3) + (-t12 * t152 + t13 * t153 - t141 * t46 - t142 * t44) * mrSges(6,2) + m(4) * (-t124 * t203 + t125 * t202 + (-t113 * t178 + t114 * t179) * qJ(3)) + ((-t178 * (Ifges(4,2) * t179 + t216) / 0.2e1 + (Ifges(4,1) * t178 + t215) * t225 + Ifges(3,5) + (-m(4) * pkin(2) - mrSges(4,1) * t179 + mrSges(4,2) * t178 - mrSges(3,1)) * pkin(7)) * qJD(2) + t221 / 0.2e1 - t234 / 0.2e1) * t185 - (t99 / 0.2e1 + t100 / 0.2e1) * t132 + (t97 / 0.2e1 - t98 / 0.2e1) * t131 + (t145 * mrSges(5,2) + t33 / 0.2e1 + t34 / 0.2e1) * t153 - (t70 / 0.2e1 + t71 / 0.2e1) * t141 + (t68 / 0.2e1 - t69 / 0.2e1) * t142 + (t31 / 0.2e1 - t32 / 0.2e1 + t145 * mrSges(5,1)) * t152; 0.2e1 * t101 * t95 + t102 * t18 + t103 * t19 + 0.2e1 * t57 * t108 + 0.2e1 * t60 * t17 + 0.2e1 * t170 * t96 + t42 * t50 + t43 * t49 + 0.2e1 * t53 * t48 + (t99 + t100) * t153 + (t97 - t98) * t152 + (t109 - t110) * t142 - (t111 + t112) * t141 + (t36 * t7 + t37 * t6 + t53 * t60) * t228 + (t101 * t57 + t197) * t229 + t197 * t230 + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * (-t116 * t141 - t117 * t142 - t152 * t66 + t153 * t67) + 0.2e1 * (t102 * t6 - t103 * t7 - t36 * t42 + t37 * t43) * mrSges(7,3) + (qJ(3) * t231 + 0.2e1 * mrSges(4,3)) * qJD(3) * (t178 ^ 2 + t179 ^ 2); m(4) * pkin(7) * t204 + m(5) * t145 + m(6) * t26 + m(7) * t14 + t133 + t38 + t39 - t5; m(6) * t57 - m(7) * t53 - t17 + t95 + t96; 0; -t232 + m(7) * (t1 * t158 + t10 * t130 + t11 * t129 + t157 * t2) + m(6) * (-pkin(4) * t13 + qJ(5) * t12 + qJD(5) * t44) + t188 + t157 * t20 + t158 * t21 + t129 * t58 + t130 * t59 + qJD(5) * t118 + qJ(5) * t61 - pkin(4) * t63 - t15 * mrSges(5,2) + t16 * mrSges(5,1) + t12 * mrSges(6,3) - t13 * mrSges(6,1); (-mrSges(5,1) - mrSges(6,1)) * t67 + (-mrSges(5,2) + mrSges(6,3)) * t66 + m(7) * (t129 * t37 + t130 * t36 + t157 * t7 + t158 * t6) + m(6) * (-pkin(4) * t67 + qJ(5) * t66 + qJD(5) * t117) + (pkin(4) * t141 - qJ(5) * t142 - qJD(5) * t152) * mrSges(6,2) + (t102 * t129 - t103 * t130 - t157 * t42 + t158 * t43) * mrSges(7,3) - t189 + t234; 0; 0.2e1 * t214 - 0.2e1 * t213 + (t129 * t158 + t130 * t157) * t228 + 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); t180 * t21 + t183 * t20 + (-t180 * t59 + t183 * t58) * qJD(6) + m(7) * (t1 * t180 + t183 * t2 + (-t10 * t180 + t11 * t183) * qJD(6)) + m(6) * t13 + t63; -t141 * mrSges(6,2) + m(7) * (t180 * t6 + t183 * t7 + (-t180 * t36 + t183 * t37) * qJD(6)) + m(6) * t67 + (t180 * t43 - t183 * t42 + (t102 * t183 + t103 * t180) * qJD(6)) * mrSges(7,3); 0; m(7) * (t129 * t180 + t130 * t183) + (m(7) * (-t157 * t180 + t158 * t183) - t233) * qJD(6); 0; -t188; t189; 0; t213 - t214; t233 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t22(1) t22(2) t22(4) t22(7) t22(11) t22(16); t22(2) t22(3) t22(5) t22(8) t22(12) t22(17); t22(4) t22(5) t22(6) t22(9) t22(13) t22(18); t22(7) t22(8) t22(9) t22(10) t22(14) t22(19); t22(11) t22(12) t22(13) t22(14) t22(15) t22(20); t22(16) t22(17) t22(18) t22(19) t22(20) t22(21);];
Mq  = res;
