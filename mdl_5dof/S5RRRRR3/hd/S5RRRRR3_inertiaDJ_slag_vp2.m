% Calculate time derivative of joint inertia matrix for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_inertiaDJ_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:11
% EndTime: 2019-07-18 17:17:17
% DurationCPUTime: 2.11s
% Computational Cost: add. (3168->291), mult. (7753->443), div. (0->0), fcn. (7347->8), ass. (0->138)
t120 = sin(qJ(4));
t117 = t120 ^ 2;
t124 = cos(qJ(4));
t161 = t124 ^ 2 + t117;
t121 = sin(qJ(3));
t122 = sin(qJ(2));
t125 = cos(qJ(3));
t126 = cos(qJ(2));
t103 = t121 * t126 + t122 * t125;
t158 = qJD(4) * t124;
t101 = t121 * t122 - t125 * t126;
t201 = qJD(2) + qJD(3);
t80 = t201 * t101;
t174 = t120 * t80;
t134 = t103 * t158 - t174;
t123 = cos(qJ(5));
t119 = sin(qJ(5));
t164 = t119 * t120;
t138 = -t123 * t124 + t164;
t200 = qJD(4) + qJD(5);
t208 = t200 * t138;
t163 = t120 * t123;
t102 = t119 * t124 + t163;
t79 = t200 * t102;
t182 = -Ifges(6,5) * t208 - Ifges(6,6) * t79;
t207 = Ifges(5,5) * t158 + t182;
t171 = t124 * t80;
t81 = t201 * t103;
t205 = -Ifges(5,5) * t171 + Ifges(5,3) * t81;
t113 = pkin(1) * t121 + pkin(5);
t177 = pkin(1) * qJD(3);
t149 = t125 * t177;
t41 = -t113 * t79 - t138 * t149;
t42 = -t102 * t149 + t113 * t208;
t204 = mrSges(6,1) * t42 - t41 * mrSges(6,2);
t67 = t79 * pkin(5);
t68 = t208 * pkin(5);
t203 = mrSges(6,1) * t68 + mrSges(6,2) * t67;
t202 = t161 * t125;
t109 = -mrSges(5,1) * t124 + mrSges(5,2) * t120;
t199 = (mrSges(5,3) * t161 - mrSges(4,2)) * t125 + (-mrSges(4,1) + t109) * t121;
t198 = 2 * m(6);
t197 = 2 * pkin(3);
t196 = -2 * Ifges(4,4);
t160 = qJD(2) * t122;
t43 = pkin(1) * t160 + pkin(2) * t81 + pkin(5) * t80;
t195 = 0.2e1 * t43;
t44 = mrSges(6,1) * t79 - mrSges(6,2) * t208;
t194 = 0.2e1 * t44;
t84 = mrSges(6,1) * t138 + mrSges(6,2) * t102;
t193 = 0.2e1 * t84;
t188 = Ifges(5,5) * t81;
t187 = pkin(1) * t125;
t184 = t79 * mrSges(6,3);
t183 = t81 * Ifges(5,6);
t181 = Ifges(5,4) * t120;
t180 = Ifges(5,4) * t124;
t179 = Ifges(5,6) * t101;
t178 = Ifges(5,6) * t120;
t176 = pkin(3) * qJD(5);
t175 = t101 * Ifges(5,5);
t172 = t124 * t43;
t169 = t103 * t120;
t168 = t103 * t124;
t167 = t103 * (pkin(3) ^ 2);
t159 = qJD(4) * t120;
t157 = qJD(5) * t119;
t156 = qJD(5) * t123;
t155 = 0.2e1 * mrSges(6,3);
t69 = -pkin(1) * t126 + pkin(2) * t101 - pkin(5) * t103;
t154 = 0.2e1 * t69;
t153 = 0.2e1 * t126;
t16 = -t103 * t79 + t138 * t80;
t17 = t102 * t80 + t103 * t208;
t152 = Ifges(6,5) * t16 + Ifges(6,6) * t17 + Ifges(6,3) * t81;
t151 = mrSges(6,3) * t176;
t150 = t123 * t208 * mrSges(6,3);
t148 = pkin(3) * t159;
t147 = t117 * t167;
t115 = -pkin(3) * t124 - pkin(2);
t145 = -t159 / 0.2e1;
t143 = pkin(3) * t81 + t172 + (-qJD(5) * t120 - t159) * t69;
t142 = mrSges(5,1) * t120 + mrSges(5,2) * t124;
t141 = Ifges(5,1) * t124 - t181;
t140 = -Ifges(5,2) * t120 + t180;
t139 = Ifges(5,5) * t120 + Ifges(5,6) * t124;
t137 = -t123 * t138 * t151 + (-pkin(3) * t184 + t102 * t151) * t119 + t207;
t62 = t102 * t103;
t63 = t138 * t103;
t36 = mrSges(6,1) * t62 - mrSges(6,2) * t63;
t56 = t103 * t140 + t179;
t136 = t197 * t36 - t179 - t56;
t55 = pkin(3) * t101 + t124 * t69;
t132 = qJD(5) * t55 + t120 * t43 + t158 * t69;
t3 = t119 * t143 + t123 * t132;
t4 = -t119 * t132 + t123 * t143;
t135 = mrSges(6,1) * t4 - t3 * mrSges(6,2) + t152;
t133 = t103 * t159 + t171;
t105 = t140 * qJD(4);
t106 = t141 * qJD(4);
t110 = Ifges(5,2) * t124 + t181;
t111 = Ifges(5,1) * t120 + t180;
t45 = -Ifges(6,4) * t208 - Ifges(6,2) * t79;
t46 = -Ifges(6,1) * t208 - Ifges(6,4) * t79;
t85 = Ifges(6,4) * t102 - Ifges(6,2) * t138;
t86 = Ifges(6,1) * t102 - Ifges(6,4) * t138;
t131 = t102 * t46 + t105 * t124 + t106 * t120 - t110 * t159 + t111 * t158 - t138 * t45 - t208 * t86 - t79 * t85;
t33 = mrSges(5,1) * t81 + mrSges(5,3) * t133;
t34 = -mrSges(5,2) * t81 - mrSges(5,3) * t134;
t70 = -mrSges(5,2) * t101 - mrSges(5,3) * t169;
t71 = mrSges(5,1) * t101 - mrSges(5,3) * t168;
t130 = -t120 * t33 + t124 * t34 + (-t120 * t70 - t124 * t71) * qJD(4);
t21 = -Ifges(5,4) * t133 - Ifges(5,2) * t134 + t183;
t22 = -Ifges(5,1) * t133 - Ifges(5,4) * t134 + t188;
t29 = t123 * t55 - t164 * t69;
t30 = t119 * t55 + t163 * t69;
t31 = -Ifges(6,4) * t63 - Ifges(6,2) * t62 + Ifges(6,6) * t101;
t32 = -Ifges(6,1) * t63 - Ifges(6,4) * t62 + Ifges(6,5) * t101;
t57 = t103 * t141 + t175;
t7 = Ifges(6,4) * t16 + Ifges(6,2) * t17 + t81 * Ifges(6,6);
t8 = Ifges(6,1) * t16 + Ifges(6,4) * t17 + t81 * Ifges(6,5);
t128 = t106 * t168 / 0.2e1 + t57 * t158 / 0.2e1 + t56 * t145 - t62 * t45 / 0.2e1 - t30 * t184 - t63 * t46 / 0.2e1 - t208 * t32 / 0.2e1 - t79 * t31 / 0.2e1 - Ifges(4,5) * t80 - Ifges(4,6) * t81 + t17 * t85 / 0.2e1 + t16 * t86 / 0.2e1 - t138 * t7 / 0.2e1 + t102 * t8 / 0.2e1 + t120 * t22 / 0.2e1 + t124 * t21 / 0.2e1 + (Ifges(6,5) * t102 - Ifges(6,6) * t138 + t139) * t81 / 0.2e1 + (-Ifges(5,6) * t159 + t207) * t101 / 0.2e1 + (-t105 / 0.2e1 + pkin(3) * t44) * t169 + (-t171 / 0.2e1 + t103 * t145) * t111 + (-t102 * t4 - t138 * t3 + t208 * t29) * mrSges(6,3) + t134 * (-t110 / 0.2e1 + pkin(3) * t84);
t114 = -pkin(2) - t187;
t108 = t115 - t187;
t107 = t121 * t177 + t148;
t104 = t142 * qJD(4);
t99 = (-mrSges(6,1) * t119 - mrSges(6,2) * t123) * t176;
t97 = t138 * pkin(5);
t96 = t102 * pkin(5);
t91 = t138 * t113;
t90 = t102 * t113;
t54 = mrSges(6,1) * t101 + mrSges(6,3) * t63;
t53 = -mrSges(6,2) * t101 - mrSges(6,3) * t62;
t26 = mrSges(5,1) * t134 - mrSges(5,2) * t133;
t11 = -mrSges(6,2) * t81 + mrSges(6,3) * t17;
t10 = mrSges(6,1) * t81 - mrSges(6,3) * t16;
t9 = -mrSges(6,1) * t17 + mrSges(6,2) * t16;
t1 = [0.2e1 * t29 * t10 + 0.2e1 * t30 * t11 + t16 * t32 + t17 * t31 + 0.2e1 * t3 * t53 + 0.2e1 * t4 * t54 - t62 * t7 - t63 * t8 + (-Ifges(6,5) * t63 - Ifges(6,6) * t62 + t103 * t196) * t81 + (-pkin(1) * (mrSges(4,1) * t81 - mrSges(4,2) * t80) + Ifges(3,4) * qJD(2) * t126) * t153 + (-t147 * t80 + t29 * t4 + t3 * t30) * t198 + m(5) * t161 * t69 * t195 - 0.2e1 * t103 * t80 * Ifges(4,1) + (-t80 * t196 + ((2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3)) * t81 + t152 + t205) * t101 + (t71 * t195 - t80 * t57 + (qJD(4) * t70 + t33) * t154 + (qJD(4) * t136 + t188 + t22) * t103) * t124 + (t70 * t195 + (-qJD(4) * t71 + t34) * t154 - t136 * t80 + (-t183 + t9 * t197 - t21 + (t124 * t167 * t198 - t175 - t57) * qJD(4)) * t103) * t120 + (0.2e1 * pkin(1) * (mrSges(4,1) * t101 + mrSges(4,2) * t103) - 0.2e1 * Ifges(3,4) * t122 + (-m(4) * pkin(1) ^ 2 + Ifges(3,1) - Ifges(3,2)) * t153) * t160; t130 * t113 + ((-t121 * t81 + t125 * t80) * mrSges(4,3) + ((-mrSges(4,3) * t101 - t120 * t71 + t124 * t70) * t125 + (mrSges(4,3) + t142) * t121 * t103) * qJD(3)) * pkin(1) + t128 + (Ifges(3,5) * t126 - Ifges(3,6) * t122) * qJD(2) + m(6) * (t29 * t42 - t3 * t91 + t30 * t41 - t4 * t90 + (-t108 * t174 + (t107 * t120 + t108 * t158) * t103) * pkin(3)) + t41 * t53 + t42 * t54 - t90 * t10 - t91 * t11 + t107 * t36 + t108 * t9 + t114 * t26; t107 * t193 + t108 * t194 + (t107 * t108 - t41 * t91 - t42 * t90) * t198 + 0.2e1 * t114 * t104 + (-t42 * t102 - t138 * t41 - t208 * t90 + t91 * t79) * t155 + 0.2e1 * (m(5) * (t113 * t202 + t114 * t121) + t199) * t177 + t131; t128 + t130 * pkin(5) + m(6) * (pkin(3) * t115 * t134 + qJD(4) * t147 + t29 * t68 - t3 * t97 - t30 * t67 - t4 * t96) - pkin(2) * t26 + t36 * t148 - t67 * t53 + t68 * t54 - t96 * t10 - t97 * t11 + t115 * t9; m(6) * (t107 * t115 + t108 * t148 - t41 * t97 - t42 * t96 + t67 * t91 - t68 * t90) + (t107 + t148) * t84 + (t108 + t115) * t44 + (-pkin(2) + t114) * t104 + (m(5) * (-pkin(2) * t121 + pkin(5) * t202) + t199) * t177 + (-(-t91 - t97) * t79 - (t90 + t96) * t208 + (-t42 - t68) * t102 - (t41 - t67) * t138) * mrSges(6,3) + t131; t148 * t193 + t115 * t194 - 0.2e1 * pkin(2) * t104 + (t115 * t148 + t67 * t97 - t68 * t96) * t198 + (-t68 * t102 + t138 * t67 - t208 * t96 + t97 * t79) * t155 + t131; mrSges(5,1) * t172 + (-mrSges(5,2) * t43 + Ifges(5,6) * t80) * t120 + (-t103 * t139 - t142 * t69) * qJD(4) + (m(6) * (t119 * t3 + t123 * t4 + t156 * t30 - t157 * t29) + t53 * t156 + t119 * t11 - t54 * t157 + t123 * t10) * pkin(3) + t135 + t205; -t142 * t149 + (t109 * t113 - t178) * qJD(4) + (t150 + m(6) * (t119 * t41 + t123 * t42 - t156 * t91 + t157 * t90)) * pkin(3) + t137 + t204; (pkin(5) * t109 - t178) * qJD(4) + (t150 + m(6) * (-t119 * t67 + t123 * t68 - t156 * t97 + t157 * t96)) * pkin(3) + t137 + t203; 0.2e1 * t99; t135; t182 + t204; t182 + t203; t99; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq  = res;
