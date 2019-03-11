% Calculate time derivative of joint inertia matrix for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:18
% EndTime: 2019-03-08 20:09:23
% DurationCPUTime: 2.20s
% Computational Cost: add. (2447->323), mult. (6187->463), div. (0->0), fcn. (5965->10), ass. (0->139)
t166 = Ifges(7,4) + Ifges(6,5);
t183 = Ifges(7,2) + Ifges(6,3);
t106 = cos(qJ(5));
t136 = qJD(5) * t106;
t103 = sin(qJ(5));
t101 = cos(pkin(11));
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t99 = sin(pkin(11));
t115 = t107 * t101 - t104 * t99;
t68 = t115 * qJD(4);
t148 = t103 * t68;
t72 = t101 * t104 + t107 * t99;
t113 = t72 * t136 + t148;
t137 = qJD(5) * t103;
t144 = t106 * t68;
t112 = t72 * t137 - t144;
t182 = m(6) + m(7);
t181 = -mrSges(6,1) - mrSges(7,1);
t180 = Ifges(7,6) * t137 + t136 * t166;
t165 = pkin(8) + qJ(3);
t82 = t165 * t99;
t83 = t165 * t101;
t179 = -t104 * t83 - t107 * t82;
t102 = cos(pkin(6));
t100 = sin(pkin(6));
t105 = sin(qJ(2));
t141 = t100 * t105;
t66 = t101 * t102 - t141 * t99;
t67 = t101 * t141 + t102 * t99;
t178 = t101 * t67 - t66 * t99;
t91 = -pkin(3) * t101 - pkin(2);
t48 = -pkin(4) * t115 - pkin(9) * t72 + t91;
t55 = -t104 * t82 + t107 * t83;
t158 = t103 * t48 + t106 * t55;
t177 = qJD(5) * t158;
t176 = m(7) * qJ(6) + mrSges(7,3);
t37 = t115 * qJD(3) + t179 * qJD(4);
t69 = t72 * qJD(4);
t47 = pkin(4) * t69 - pkin(9) * t68;
t5 = -t103 * t37 + t106 * t47 - t177;
t175 = -2 * mrSges(5,3);
t174 = -2 * Ifges(5,4);
t38 = qJD(3) * t72 + qJD(4) * t55;
t173 = 0.2e1 * t38;
t172 = -0.2e1 * t179;
t170 = t38 * t179;
t108 = cos(qJ(2));
t138 = qJD(2) * t108;
t128 = t100 * t138;
t40 = t104 * t66 + t107 * t67;
t23 = qJD(4) * t40 + t128 * t72;
t39 = t104 * t67 - t107 * t66;
t11 = t39 * t23;
t168 = t115 * Ifges(6,6);
t139 = qJD(2) * t100;
t129 = t105 * t139;
t22 = -qJD(4) * t39 + t115 * t128;
t140 = t100 * t108;
t30 = t103 * t40 + t106 * t140;
t9 = -qJD(5) * t30 + t103 * t129 + t106 * t22;
t167 = t9 * t106;
t26 = mrSges(6,1) * t69 + mrSges(6,3) * t112;
t27 = -t69 * mrSges(7,1) - mrSges(7,2) * t112;
t164 = t27 - t26;
t28 = -mrSges(6,2) * t69 - mrSges(6,3) * t113;
t29 = -mrSges(7,2) * t113 + mrSges(7,3) * t69;
t163 = t28 + t29;
t151 = Ifges(7,5) * t106;
t119 = Ifges(7,3) * t103 + t151;
t32 = -Ifges(7,6) * t115 + t119 * t72;
t153 = Ifges(6,4) * t106;
t120 = -Ifges(6,2) * t103 + t153;
t33 = t120 * t72 - t168;
t162 = t32 - t33;
t152 = Ifges(7,5) * t103;
t121 = Ifges(7,1) * t106 + t152;
t34 = -Ifges(7,4) * t115 + t121 * t72;
t154 = Ifges(6,4) * t103;
t122 = Ifges(6,1) * t106 - t154;
t35 = -Ifges(6,5) * t115 + t122 * t72;
t161 = t34 + t35;
t147 = t103 * t72;
t49 = mrSges(6,2) * t115 - mrSges(6,3) * t147;
t52 = -mrSges(7,2) * t147 - mrSges(7,3) * t115;
t160 = t49 + t52;
t143 = t106 * t72;
t50 = -mrSges(6,1) * t115 - mrSges(6,3) * t143;
t51 = mrSges(7,1) * t115 + mrSges(7,2) * t143;
t159 = -t50 + t51;
t85 = -t106 * mrSges(6,1) + t103 * mrSges(6,2);
t157 = t85 - mrSges(5,1);
t155 = t101 ^ 2 + t99 ^ 2;
t130 = t103 * t140;
t10 = -qJD(5) * t130 + t103 * t22 - t106 * t129 + t136 * t40;
t150 = t10 * t103;
t145 = t105 * t100 ^ 2;
t142 = qJD(5) * t72;
t135 = qJD(6) * t106;
t86 = -Ifges(7,3) * t106 + t152;
t87 = Ifges(6,2) * t106 + t154;
t134 = t86 / 0.2e1 - t87 / 0.2e1;
t88 = Ifges(7,1) * t103 - t151;
t89 = Ifges(6,1) * t103 + t153;
t133 = t88 / 0.2e1 + t89 / 0.2e1;
t44 = t69 * mrSges(5,1) + t68 * mrSges(5,2);
t125 = -t179 * t23 + t38 * t39;
t124 = t103 * mrSges(6,1) + t106 * mrSges(6,2);
t84 = -t106 * mrSges(7,1) - t103 * mrSges(7,3);
t123 = t103 * mrSges(7,1) - t106 * mrSges(7,3);
t118 = -pkin(5) * t106 - qJ(6) * t103;
t20 = -t103 * t55 + t106 * t48;
t114 = t113 * Ifges(7,6) + t144 * t166 + t183 * t69;
t4 = t103 * t47 + t106 * t37 + t48 * t136 - t137 * t55;
t79 = -pkin(4) + t118;
t78 = t122 * qJD(5);
t77 = t121 * qJD(5);
t76 = t120 * qJD(5);
t75 = t119 * qJD(5);
t74 = t124 * qJD(5);
t73 = t123 * qJD(5);
t65 = -pkin(5) * t137 + qJ(6) * t136 + qJD(6) * t103;
t46 = t124 * t72;
t45 = t123 * t72;
t31 = t106 * t40 - t130;
t24 = (pkin(5) * t103 - qJ(6) * t106) * t72 - t179;
t19 = mrSges(6,1) * t113 - mrSges(6,2) * t112;
t18 = mrSges(7,1) * t113 + mrSges(7,3) * t112;
t17 = pkin(5) * t115 - t20;
t16 = -qJ(6) * t115 + t158;
t15 = -Ifges(6,1) * t112 - Ifges(6,4) * t113 + t69 * Ifges(6,5);
t14 = -Ifges(7,1) * t112 + t69 * Ifges(7,4) + Ifges(7,5) * t113;
t13 = -Ifges(6,4) * t112 - Ifges(6,2) * t113 + t69 * Ifges(6,6);
t12 = -Ifges(7,5) * t112 + t69 * Ifges(7,6) + Ifges(7,3) * t113;
t8 = (pkin(5) * t68 + qJ(6) * t142) * t103 + (-qJ(6) * t68 + (pkin(5) * qJD(5) - qJD(6)) * t72) * t106 + t38;
t6 = pkin(9) * t167;
t3 = -pkin(5) * t69 - t5;
t1 = qJ(6) * t69 - qJD(6) * t115 + t4;
t2 = [0.2e1 * m(4) * (t178 * t100 - t145) * t138 + 0.2e1 * t182 * (t10 * t30 + t31 * t9 + t11) + 0.2e1 * (-t138 * t145 + t40 * t22 + t11) * m(5); t160 * t9 + (t18 + t19) * t39 + t163 * t31 + t164 * t30 + (t45 + t46) * t23 + t159 * t10 + (t115 * t22 + t23 * t72 + t39 * t68 - t40 * t69) * mrSges(5,3) + (-t108 * t44 + ((mrSges(4,3) * t155 - mrSges(3,2)) * t108 + (-mrSges(4,1) * t101 - mrSges(5,1) * t115 + mrSges(4,2) * t99 + mrSges(5,2) * t72 - mrSges(3,1)) * t105) * qJD(2)) * t100 + m(7) * (t1 * t31 + t10 * t17 + t16 * t9 + t23 * t24 + t3 * t30 + t39 * t8) + m(6) * (-t10 * t20 + t158 * t9 - t30 * t5 + t31 * t4 + t125) + m(5) * (t129 * t91 + t22 * t55 + t37 * t40 + t125) + m(4) * (t178 * qJD(3) + (qJ(3) * t108 * t155 - pkin(2) * t105) * t139); t55 * t69 * t175 + 0.2e1 * t1 * t52 + 0.2e1 * t16 * t29 + 0.2e1 * t17 * t27 + 0.2e1 * t24 * t18 + t19 * t172 + 0.2e1 * t20 * t26 + 0.2e1 * t158 * t28 + 0.2e1 * t3 * t51 + t46 * t173 + 0.2e1 * t4 * t49 + 0.2e1 * t91 * t44 + 0.2e1 * t8 * t45 + 0.2e1 * t5 * t50 + (mrSges(5,3) * t172 + t103 * t162 + t106 * t161) * t68 + 0.2e1 * m(7) * (t1 * t16 + t17 * t3 + t24 * t8) + 0.2e1 * m(6) * (t158 * t4 + t20 * t5 - t170) + 0.2e1 * m(5) * (t37 * t55 - t170) - (t37 * t175 + (-Ifges(6,6) * t103 + t174) * t68 + ((2 * Ifges(5,2)) + t183) * t69 + t114) * t115 + (mrSges(5,3) * t173 + 0.2e1 * Ifges(5,1) * t68 + (t14 + t15) * t106 + (t12 - t13) * t103 + (t174 + t166 * t106 + (-Ifges(6,6) + Ifges(7,6)) * t103) * t69 + ((t162 + t168) * t106 + (t115 * t166 - t161) * t103) * qJD(5)) * t72 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3) * t155; (m(4) + m(5)) * t129 + t182 * (-t10 * t106 + t103 * t9 + t31 * t136 + t137 * t30); -t164 * t106 + t163 * t103 + (t159 * t103 + t160 * t106) * qJD(5) + m(7) * (t1 * t103 - t106 * t3 + (t103 * t17 + t106 * t16) * qJD(5)) + m(6) * (t103 * t4 + t106 * t5 + (-t103 * t20 + t106 * t158) * qJD(5)) + t44; 0; -t22 * mrSges(5,2) + (t73 + t74) * t39 + (t84 + t157) * t23 + m(6) * (-pkin(4) * t23 + t6) + m(7) * (t23 * t79 - t39 * t65 + t6) + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t136 * t30 - t137 * t31 + t150) * pkin(9) + (mrSges(6,3) + mrSges(7,2)) * (t150 + t167 + (-t103 * t31 + t106 * t30) * qJD(5)); -t37 * mrSges(5,2) - pkin(4) * t19 + m(7) * (-t24 * t65 + t79 * t8) - t65 * t45 + Ifges(5,5) * t68 - Ifges(5,6) * t69 + t24 * t73 - t179 * t74 + t79 * t18 + t8 * t84 + (-m(6) * pkin(4) + t157) * t38 + (t3 * mrSges(7,2) - t5 * mrSges(6,3) + t14 / 0.2e1 + t15 / 0.2e1 + (t75 / 0.2e1 - t76 / 0.2e1) * t72 + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t69 + t134 * t68 + (-t16 * mrSges(7,2) - t158 * mrSges(6,3) + t32 / 0.2e1 - t33 / 0.2e1 + t168 / 0.2e1 - t133 * t72) * qJD(5) + (-t160 * qJD(5) + m(6) * (-t5 - t177) + m(7) * (-t16 * qJD(5) + t3) + t164) * pkin(9)) * t103 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) - t12 / 0.2e1 + t13 / 0.2e1 + (t77 / 0.2e1 + t78 / 0.2e1) * t72 + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t69 + t133 * t68 + (t17 * mrSges(7,2) - t20 * mrSges(6,3) + t34 / 0.2e1 + t35 / 0.2e1 + t134 * t72) * qJD(5) + (t159 * qJD(5) + m(7) * (t17 * qJD(5) + t1) + m(6) * (-t20 * qJD(5) + t4) + t163) * pkin(9)) * t106 - t180 * t115 / 0.2e1; 0; -0.2e1 * pkin(4) * t74 + 0.2e1 * t73 * t79 + 0.2e1 * (-m(7) * t79 - t84) * t65 + (-t75 + t76) * t106 + (t77 + t78) * t103 + ((t88 + t89) * t106 + (t86 - t87) * t103) * qJD(5); m(7) * qJD(6) * t31 + (-mrSges(6,2) + t176) * t9 + (-m(7) * pkin(5) + t181) * t10; -Ifges(6,6) * t148 - pkin(5) * t27 + m(7) * (-pkin(5) * t3 + qJ(6) * t1 + qJD(6) * t16) + qJD(6) * t52 + qJ(6) * t29 + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t5 * mrSges(6,1) - t3 * mrSges(7,1) + (-Ifges(6,6) * t106 - t103 * t166) * t142 + t114; m(7) * t65 + ((-mrSges(6,2) + mrSges(7,3)) * t106 + t181 * t103) * qJD(5); -Ifges(6,6) * t137 + (qJD(5) * t118 + t135) * mrSges(7,2) + (m(7) * t135 + (m(7) * t118 + t84 + t85) * qJD(5)) * pkin(9) + t180; 0.2e1 * t176 * qJD(6); m(7) * t10; m(7) * t3 + t27; m(7) * t137; (m(7) * pkin(9) + mrSges(7,2)) * t136; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
