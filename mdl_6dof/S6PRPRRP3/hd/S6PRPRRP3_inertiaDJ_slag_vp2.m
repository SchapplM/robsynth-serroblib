% Calculate time derivative of joint inertia matrix for
% S6PRPRRP3
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:04:46
% EndTime: 2019-03-08 20:04:50
% DurationCPUTime: 1.96s
% Computational Cost: add. (2432->326), mult. (6168->470), div. (0->0), fcn. (5931->10), ass. (0->148)
t167 = Ifges(6,5) + Ifges(7,5);
t185 = Ifges(6,3) + Ifges(7,3);
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t138 = qJD(5) * t106;
t101 = cos(pkin(11));
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t99 = sin(pkin(11));
t72 = t101 * t104 + t107 * t99;
t129 = t72 * t138;
t113 = t107 * t101 - t104 * t99;
t68 = t113 * qJD(4);
t112 = t103 * t68 + t129;
t139 = qJD(5) * t103;
t145 = t106 * t68;
t111 = t72 * t139 - t145;
t184 = (-m(6) * pkin(9) - mrSges(6,3));
t165 = pkin(8) + qJ(3);
t79 = t165 * t99;
t80 = t165 * t101;
t183 = -t104 * t80 - t107 * t79;
t102 = cos(pkin(6));
t100 = sin(pkin(6));
t105 = sin(qJ(2));
t143 = t100 * t105;
t66 = t101 * t102 - t143 * t99;
t67 = t101 * t143 + t102 * t99;
t182 = t101 * t67 - t66 * t99;
t181 = -m(6) * pkin(4) - mrSges(6,1) * t106 + mrSges(6,2) * t103 - mrSges(5,1);
t82 = -mrSges(7,1) * t106 + mrSges(7,2) * t103;
t91 = -pkin(5) * t106 - pkin(4);
t180 = m(7) * t91 + t82;
t179 = 2 * m(5);
t178 = 0.2e1 * m(7);
t177 = -2 * mrSges(5,3);
t176 = -2 * mrSges(7,3);
t54 = -t104 * t79 + t107 * t80;
t34 = qJD(3) * t72 + qJD(4) * t54;
t175 = 0.2e1 * t34;
t174 = -0.2e1 * t183;
t171 = m(7) * pkin(5);
t170 = t34 * t183;
t108 = cos(qJ(2));
t140 = qJD(2) * t108;
t127 = t100 * t140;
t37 = t104 * t66 + t107 * t67;
t21 = qJD(4) * t37 + t127 * t72;
t36 = t104 * t67 - t107 * t66;
t169 = t36 * t21;
t166 = -Ifges(6,6) - Ifges(7,6);
t164 = -qJ(6) - pkin(9);
t69 = t72 * qJD(4);
t22 = mrSges(7,1) * t69 + mrSges(7,3) * t111;
t23 = mrSges(6,1) * t69 + mrSges(6,3) * t111;
t163 = t22 + t23;
t24 = -mrSges(7,2) * t69 - mrSges(7,3) * t112;
t25 = -mrSges(6,2) * t69 - mrSges(6,3) * t112;
t162 = t24 + t25;
t151 = Ifges(7,4) * t106;
t115 = -Ifges(7,2) * t103 + t151;
t28 = -Ifges(7,6) * t113 + t115 * t72;
t153 = Ifges(6,4) * t106;
t116 = -Ifges(6,2) * t103 + t153;
t29 = -Ifges(6,6) * t113 + t116 * t72;
t161 = -t28 - t29;
t152 = Ifges(7,4) * t103;
t117 = Ifges(7,1) * t106 - t152;
t30 = -Ifges(7,5) * t113 + t117 * t72;
t154 = Ifges(6,4) * t103;
t118 = Ifges(6,1) * t106 - t154;
t31 = -Ifges(6,5) * t113 + t118 * t72;
t160 = t30 + t31;
t148 = t103 * t72;
t48 = mrSges(7,2) * t113 - mrSges(7,3) * t148;
t49 = mrSges(6,2) * t113 - mrSges(6,3) * t148;
t159 = t48 + t49;
t144 = t106 * t72;
t50 = -mrSges(7,1) * t113 - mrSges(7,3) * t144;
t51 = -mrSges(6,1) * t113 - mrSges(6,3) * t144;
t158 = t50 + t51;
t90 = -pkin(3) * t101 - pkin(2);
t47 = -pkin(4) * t113 - pkin(9) * t72 + t90;
t52 = t106 * t54;
t19 = t103 * t47 + t52;
t73 = mrSges(7,1) * t139 + mrSges(7,2) * t138;
t156 = t101 ^ 2 + t99 ^ 2;
t155 = mrSges(6,2) * t106;
t146 = t105 * t100 ^ 2;
t142 = t100 * t108;
t141 = qJD(2) * t100;
t137 = m(6) / 0.2e1 + m(7) / 0.2e1;
t33 = t113 * qJD(3) + t183 * qJD(4);
t46 = pkin(4) * t69 - pkin(9) * t68;
t136 = t103 * t46 + t106 * t33 + t47 * t138;
t135 = pkin(5) * t139;
t134 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t85 = Ifges(7,2) * t106 + t152;
t86 = Ifges(6,2) * t106 + t154;
t133 = -t85 / 0.2e1 - t86 / 0.2e1;
t87 = Ifges(7,1) * t103 + t151;
t88 = Ifges(6,1) * t103 + t153;
t132 = t87 / 0.2e1 + t88 / 0.2e1;
t131 = mrSges(7,1) + t171;
t128 = t105 * t141;
t43 = t69 * mrSges(5,1) + t68 * mrSges(5,2);
t126 = t166 * t103;
t125 = -t103 * t33 + t106 * t46;
t18 = -t103 * t54 + t106 * t47;
t124 = qJD(5) * t164;
t122 = t167 * t145 + t185 * t69;
t121 = -t183 * t21 + t34 * t36;
t120 = -(2 * Ifges(5,4)) + t126;
t119 = mrSges(6,1) * t103 + t155;
t114 = -qJ(6) * t68 - qJD(6) * t72;
t26 = -t103 * t37 - t106 * t142;
t110 = t103 * t142 - t106 * t37;
t16 = t112 * mrSges(7,1) - t111 * mrSges(7,2);
t95 = Ifges(6,5) * t138;
t94 = Ifges(7,5) * t138;
t84 = t164 * t106;
t81 = t164 * t103;
t78 = t118 * qJD(5);
t77 = t117 * qJD(5);
t76 = t116 * qJD(5);
t75 = t115 * qJD(5);
t74 = t119 * qJD(5);
t65 = -qJD(6) * t103 + t106 * t124;
t64 = qJD(6) * t106 + t103 * t124;
t45 = t119 * t72;
t44 = (mrSges(7,1) * t103 + mrSges(7,2) * t106) * t72;
t35 = pkin(5) * t148 - t183;
t20 = -qJD(4) * t36 + t113 * t127;
t17 = mrSges(6,1) * t112 - mrSges(6,2) * t111;
t15 = pkin(5) * t112 + t34;
t14 = -qJ(6) * t148 + t19;
t13 = -Ifges(6,1) * t111 - Ifges(6,4) * t112 + t69 * Ifges(6,5);
t12 = -Ifges(7,1) * t111 - Ifges(7,4) * t112 + t69 * Ifges(7,5);
t11 = -Ifges(6,4) * t111 - Ifges(6,2) * t112 + t69 * Ifges(6,6);
t10 = -Ifges(7,4) * t111 - Ifges(7,2) * t112 + t69 * Ifges(7,6);
t9 = -pkin(5) * t113 - qJ(6) * t144 + t18;
t8 = qJD(5) * t110 - t103 * t20 + t106 * t128;
t7 = qJD(5) * t26 + t103 * t128 + t106 * t20;
t6 = -qJD(5) * t19 + t125;
t5 = -t139 * t54 + t136;
t4 = -qJ(6) * t129 + (-qJD(5) * t54 + t114) * t103 + t136;
t3 = pkin(5) * t69 + t114 * t106 + (-t52 + (qJ(6) * t72 - t47) * t103) * qJD(5) + t125;
t1 = [(t37 * t20 + t169) * t179 + 0.4e1 * t137 * (-t110 * t7 + t26 * t8 + t169) + 0.4e1 * (-m(5) * t146 / 0.2e1 + m(4) * (t100 * t182 - t146) / 0.2e1) * t140; t158 * t8 + t159 * t7 + (t16 + t17) * t36 - t162 * t110 + t163 * t26 + (t44 + t45) * t21 + (t113 * t20 + t21 * t72 + t36 * t68 - t37 * t69) * mrSges(5,3) + (-t108 * t43 + ((mrSges(4,3) * t156 - mrSges(3,2)) * t108 + (-mrSges(4,1) * t101 - mrSges(5,1) * t113 + mrSges(4,2) * t99 + mrSges(5,2) * t72 - mrSges(3,1)) * t105) * qJD(2)) * t100 + m(7) * (-t110 * t4 + t14 * t7 + t15 * t36 + t21 * t35 + t26 * t3 + t8 * t9) + m(6) * (-t110 * t5 + t18 * t8 + t19 * t7 + t26 * t6 + t121) + m(5) * (t128 * t90 + t20 * t54 + t33 * t37 + t121) + m(4) * (t182 * qJD(3) + (qJ(3) * t108 * t156 - pkin(2) * t105) * t141); t54 * t69 * t177 + 0.2e1 * t14 * t24 + 0.2e1 * t15 * t44 + 0.2e1 * t35 * t16 + t17 * t174 + 0.2e1 * t18 * t23 + 0.2e1 * t19 * t25 + 0.2e1 * t9 * t22 + 0.2e1 * t3 * t50 + t45 * t175 + 0.2e1 * t4 * t48 + 0.2e1 * t90 * t43 + 0.2e1 * t5 * t49 + 0.2e1 * t6 * t51 + (mrSges(5,3) * t174 + t103 * t161 + t106 * t160) * t68 + (t14 * t4 + t15 * t35 + t3 * t9) * t178 + 0.2e1 * m(6) * (t18 * t6 + t19 * t5 - t170) + (t33 * t54 - t170) * t179 - (t33 * t177 + ((2 * Ifges(5,2)) + t185) * t69 + t120 * t68 + t122) * t113 + (mrSges(5,3) * t175 + 0.2e1 * Ifges(5,1) * t68 + (t12 + t13) * t106 + (-t10 - t11) * t103 + (t106 * t167 + t120) * t69 + ((-t113 * t166 + t161) * t106 + (t113 * t167 - t160) * t103) * qJD(5)) * t72 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * t156 * qJD(3); (m(4) + m(5)) * t128 + 0.2e1 * t137 * (t103 * t7 + t106 * t8 + (-t103 * t26 - t106 * t110) * qJD(5)); t163 * t106 + t162 * t103 + (-t158 * t103 + t159 * t106) * qJD(5) + m(7) * (t103 * t4 + t106 * t3 + (-t103 * t9 + t106 * t14) * qJD(5)) + m(6) * (t103 * t5 + t106 * t6 + (-t103 * t18 + t106 * t19) * qJD(5)) + t43; 0; -t20 * mrSges(5,2) + (t73 + t74) * t36 + m(7) * (-t110 * t64 + t135 * t36 + t26 * t65 - t7 * t84 + t8 * t81) + (t180 + t181) * t21 + (mrSges(7,3) - t184) * (-t8 * t103 + t7 * t106 + (t103 * t110 - t106 * t26) * qJD(5)); m(7) * (t14 * t64 + t15 * t91 + t3 * t81 - t4 * t84 + t65 * t9) - pkin(4) * t17 - t33 * mrSges(5,2) + t64 * t48 + t65 * t50 + Ifges(5,5) * t68 - Ifges(5,6) * t69 + t35 * t73 - t183 * t74 + t81 * t22 + t15 * t82 - t84 * t24 + t91 * t16 - (t94 / 0.2e1 + t95 / 0.2e1) * t113 + t181 * t34 + (-t3 * mrSges(7,3) - t6 * mrSges(6,3) + t12 / 0.2e1 + t13 / 0.2e1 + (-t75 / 0.2e1 - t76 / 0.2e1) * t72 + (Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1) * t69 + t133 * t68 + (-m(6) * t6 - t23) * pkin(9) + (pkin(5) * t44 - pkin(9) * t49 - t14 * mrSges(7,3) - t28 / 0.2e1 - t29 / 0.2e1 - t132 * t72 - t134 * t113 + t35 * t171 + t184 * t19) * qJD(5)) * t103 + (t4 * mrSges(7,3) + t5 * mrSges(6,3) + t10 / 0.2e1 + t11 / 0.2e1 + (t77 / 0.2e1 + t78 / 0.2e1) * t72 - t134 * t69 + t132 * t68 + (m(6) * t5 + t25) * pkin(9) + (-t9 * mrSges(7,3) - t18 * mrSges(6,3) + t30 / 0.2e1 + t31 / 0.2e1 + t133 * t72 + (-m(6) * t18 - t51) * pkin(9)) * qJD(5)) * t106; m(7) * (t103 * t64 + t106 * t65 + (-t103 * t81 - t106 * t84) * qJD(5)); -0.2e1 * pkin(4) * t74 + 0.2e1 * t91 * t73 + (-t64 * t84 + t65 * t81) * t178 + (t65 * t176 + t77 + t78 + (0.2e1 * pkin(5) * t180 - t84 * t176 - t85 - t86) * qJD(5)) * t103 + (0.2e1 * t64 * mrSges(7,3) + t75 + t76 + (t176 * t81 + t87 + t88) * qJD(5)) * t106; (-mrSges(6,2) - mrSges(7,2)) * t7 + (mrSges(6,1) + t131) * t8; mrSges(6,1) * t6 + mrSges(7,1) * t3 - mrSges(6,2) * t5 - mrSges(7,2) * t4 + t68 * t126 + (m(7) * t3 + t22) * pkin(5) + (-t167 * t103 + t166 * t106) * t72 * qJD(5) + t122; (-t155 + (-mrSges(6,1) - t171) * t103) * qJD(5) - t73; -mrSges(7,2) * t64 + t94 + t95 + t131 * t65 + ((-mrSges(6,1) * pkin(9) - mrSges(7,3) * pkin(5)) * t106 + (mrSges(6,2) * pkin(9) + t166) * t103) * qJD(5); 0; m(7) * t21; m(7) * t15 + t16; 0; m(7) * t135 + t73; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
