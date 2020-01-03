% Calculate time derivative of joint inertia matrix for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR15_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR15_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:12
% EndTime: 2019-12-31 20:41:17
% DurationCPUTime: 1.79s
% Computational Cost: add. (2027->282), mult. (4358->432), div. (0->0), fcn. (3393->6), ass. (0->142)
t161 = pkin(3) + pkin(6);
t104 = cos(qJ(4));
t101 = sin(qJ(4));
t105 = cos(qJ(2));
t138 = qJD(4) * t105;
t132 = t101 * t138;
t102 = sin(qJ(2));
t142 = qJD(2) * t102;
t110 = t104 * t142 + t132;
t173 = m(4) * pkin(6);
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t125 = pkin(2) * t142 - qJD(3) * t102;
t51 = (pkin(7) * t102 - qJ(3) * t105) * qJD(2) + t125;
t141 = qJD(2) * t105;
t77 = t161 * t141;
t126 = -t101 * t51 + t104 * t77;
t106 = -pkin(2) - pkin(7);
t145 = qJ(3) * t102;
t62 = t106 * t105 - pkin(1) - t145;
t127 = pkin(8) * t105 - t62;
t84 = t161 * t102;
t69 = t101 * t84;
t7 = (-pkin(8) * t101 * t102 + pkin(4) * t105) * qJD(2) + (t127 * t104 - t69) * qJD(4) + t126;
t139 = qJD(4) * t104;
t140 = qJD(4) * t101;
t11 = t101 * t77 + t104 * t51 + t84 * t139 - t62 * t140;
t8 = t110 * pkin(8) + t11;
t70 = t104 * t84;
t27 = pkin(4) * t102 + t127 * t101 + t70;
t143 = t104 * t105;
t37 = t104 * t62 + t69;
t31 = -pkin(8) * t143 + t37;
t9 = -t100 * t31 + t103 * t27;
t2 = t9 * qJD(5) + t100 * t7 + t103 * t8;
t10 = t100 * t27 + t103 * t31;
t3 = -t10 * qJD(5) - t100 * t8 + t103 * t7;
t172 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t12 = -t37 * qJD(4) + t126;
t171 = t101 * t11 + t104 * t12;
t170 = qJD(4) + qJD(5);
t144 = t103 * t104;
t112 = t100 * t101 - t144;
t113 = t100 * t104 + t103 * t101;
t137 = qJD(5) * t100;
t38 = -t100 * t140 - t101 * t137 + t170 * t144;
t39 = t170 * t113;
t169 = (-t100 * t112 - t103 * t113) * qJD(5) - t100 * t38 + t103 * t39;
t168 = t105 * t170;
t167 = 2 * m(6);
t166 = -0.2e1 * pkin(1);
t57 = -qJ(3) * t141 + t125;
t165 = 0.2e1 * t57;
t119 = -pkin(2) * t105 - t145;
t80 = -pkin(1) + t119;
t164 = -0.2e1 * t80;
t163 = -t113 / 0.2e1;
t162 = -t112 / 0.2e1;
t159 = -t101 / 0.2e1;
t158 = t102 / 0.2e1;
t157 = -t104 / 0.2e1;
t156 = pkin(8) - t106;
t155 = -Ifges(6,5) * t39 - Ifges(6,6) * t38;
t85 = t161 * t105;
t154 = Ifges(5,4) * t101;
t153 = Ifges(5,4) * t104;
t152 = Ifges(5,5) * t102;
t122 = Ifges(5,1) * t101 + t153;
t53 = -t122 * t105 + t152;
t150 = t101 * t53;
t83 = Ifges(5,1) * t104 - t154;
t149 = t101 * t83;
t121 = Ifges(5,2) * t104 + t154;
t52 = t102 * Ifges(5,6) - t121 * t105;
t147 = t104 * t52;
t82 = -Ifges(5,2) * t101 + t153;
t146 = t104 * t82;
t136 = qJD(5) * t103;
t135 = 2 * mrSges(6,3);
t23 = t112 * t168 + t113 * t142;
t24 = -t112 * t142 + t113 * t168;
t134 = Ifges(6,5) * t23 + Ifges(6,6) * t24 + Ifges(6,3) * t141;
t133 = mrSges(4,1) + t173;
t131 = t104 * t138;
t130 = t101 * t142;
t128 = -t39 * mrSges(6,1) - t38 * mrSges(6,2);
t79 = t156 * t104;
t124 = -t112 * t39 - t113 * t38;
t123 = mrSges(5,1) * t104 - mrSges(5,2) * t101;
t81 = mrSges(5,1) * t101 + mrSges(5,2) * t104;
t120 = -Ifges(5,5) * t101 - Ifges(5,6) * t104;
t78 = t156 * t101;
t44 = -t100 * t79 - t103 * t78;
t43 = t100 * t78 - t103 * t79;
t36 = -t101 * t62 + t70;
t116 = t101 * t36 - t104 * t37;
t74 = mrSges(5,3) * t101 * t105 + mrSges(5,1) * t102;
t75 = -mrSges(5,2) * t102 - mrSges(5,3) * t143;
t115 = -t101 * t74 + t104 * t75;
t65 = t156 * t140;
t66 = qJD(4) * t79;
t16 = t43 * qJD(5) + t100 * t65 - t103 * t66;
t17 = -t44 * qJD(5) + t100 * t66 + t103 * t65;
t114 = t17 * mrSges(6,1) - t16 * mrSges(6,2) + t155;
t111 = Ifges(5,5) * t130 + t110 * Ifges(5,6) + Ifges(5,3) * t141 + t134;
t109 = t130 - t131;
t108 = t10 * t38 - t112 * t3 + t113 * t2 - t39 * t9;
t107 = t112 * t17 - t113 * t16 - t38 * t44 + t39 * t43;
t92 = pkin(4) * t101 + qJ(3);
t89 = pkin(4) * t139 + qJD(3);
t76 = t161 * t142;
t73 = t122 * qJD(4);
t72 = t121 * qJD(4);
t71 = t123 * qJD(4);
t63 = (-mrSges(6,1) * t100 - mrSges(6,2) * t103) * qJD(5) * pkin(4);
t61 = t123 * t105;
t58 = pkin(4) * t143 + t85;
t55 = t113 * t105;
t54 = t112 * t105;
t50 = mrSges(5,1) * t141 - t109 * mrSges(5,3);
t49 = -mrSges(5,2) * t141 + t110 * mrSges(5,3);
t48 = mrSges(6,1) * t102 + mrSges(6,3) * t55;
t47 = -mrSges(6,2) * t102 + mrSges(6,3) * t54;
t45 = -pkin(4) * t132 + (-pkin(4) * t104 - t161) * t142;
t42 = -Ifges(6,1) * t112 - Ifges(6,4) * t113;
t41 = -Ifges(6,4) * t112 - Ifges(6,2) * t113;
t40 = mrSges(6,1) * t113 - mrSges(6,2) * t112;
t32 = -t110 * mrSges(5,1) + t109 * mrSges(5,2);
t30 = -mrSges(6,1) * t54 - mrSges(6,2) * t55;
t29 = -t83 * t138 + (Ifges(5,5) * t105 + t122 * t102) * qJD(2);
t28 = -t82 * t138 + (Ifges(5,6) * t105 + t121 * t102) * qJD(2);
t26 = -Ifges(6,1) * t55 + Ifges(6,4) * t54 + Ifges(6,5) * t102;
t25 = -Ifges(6,4) * t55 + Ifges(6,2) * t54 + Ifges(6,6) * t102;
t20 = -Ifges(6,1) * t39 - Ifges(6,4) * t38;
t19 = -Ifges(6,4) * t39 - Ifges(6,2) * t38;
t18 = mrSges(6,1) * t38 - mrSges(6,2) * t39;
t14 = -mrSges(6,2) * t141 + mrSges(6,3) * t24;
t13 = mrSges(6,1) * t141 - mrSges(6,3) * t23;
t6 = -mrSges(6,1) * t24 + mrSges(6,2) * t23;
t5 = Ifges(6,1) * t23 + Ifges(6,4) * t24 + Ifges(6,5) * t141;
t4 = Ifges(6,4) * t23 + Ifges(6,2) * t24 + Ifges(6,6) * t141;
t1 = [0.2e1 * t85 * t32 + 0.2e1 * t12 * t74 + 0.2e1 * t11 * t75 - 0.2e1 * t76 * t61 - t55 * t5 + 0.2e1 * t58 * t6 + 0.2e1 * t45 * t30 + 0.2e1 * t2 * t47 + 0.2e1 * t3 * t48 + 0.2e1 * t37 * t49 + 0.2e1 * t36 * t50 + t54 * t4 + t24 * t25 + t23 * t26 + 0.2e1 * t9 * t13 + 0.2e1 * t10 * t14 + m(4) * t80 * t165 + 0.2e1 * m(5) * (t11 * t37 + t12 * t36 - t76 * t85) + (t10 * t2 + t3 * t9 + t45 * t58) * t167 + (-0.2e1 * t57 * mrSges(4,3) + t111) * t102 + (mrSges(4,2) * t165 - t101 * t29 - t104 * t28 + (t101 * t52 + (-t53 - t152) * t104) * qJD(4)) * t105 + ((mrSges(3,1) * t166 + mrSges(4,2) * t164 + t150 + t147 + 0.2e1 * (-Ifges(4,6) - Ifges(3,4)) * t102) * t102 + (mrSges(3,2) * t166 + mrSges(4,3) * t164 - Ifges(6,5) * t55 + Ifges(6,6) * t54 + (0.2e1 * Ifges(3,4) + 0.2e1 * Ifges(4,6) + t120) * t105 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(5,3) + Ifges(6,3)) * t102) * t105) * qJD(2); -t108 * mrSges(6,3) + (t133 * qJD(3) - t72 * t157 - t73 * t159) * t105 + (t106 * t50 - t12 * mrSges(5,3) + t29 / 0.2e1) * t104 + (t106 * t49 - t11 * mrSges(5,3) - t28 / 0.2e1) * t101 + t85 * t71 + t89 * t30 + t92 * t6 - t76 * t81 - t55 * t20 / 0.2e1 + t58 * t18 + qJD(3) * t61 + t44 * t14 + t45 * t40 + t16 * t47 + t17 * t48 + t54 * t19 / 0.2e1 - t38 * t25 / 0.2e1 - t39 * t26 / 0.2e1 + t24 * t41 / 0.2e1 + t23 * t42 / 0.2e1 + t43 * t13 + qJ(3) * t32 + m(6) * (t10 * t16 + t17 * t9 + t2 * t44 + t3 * t43 + t45 * t92 + t58 * t89) + t155 * t158 + t5 * t162 + t4 * t163 + (t119 * t173 + (Ifges(5,5) * t104 / 0.2e1 + Ifges(5,6) * t159 + Ifges(6,5) * t162 + Ifges(6,6) * t163 - Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,1) + (-mrSges(3,1) + mrSges(4,2)) * pkin(6)) * t105 + (Ifges(4,5) - Ifges(3,6) + t146 / 0.2e1 + t149 / 0.2e1 - qJ(3) * mrSges(4,1) + (mrSges(3,2) - mrSges(4,3)) * pkin(6)) * t102) * qJD(2) + m(5) * (-qJ(3) * t76 + qJD(3) * t85 + t171 * t106) + (-t147 / 0.2e1 - t150 / 0.2e1 + t120 * t158 + (t83 * t157 + t101 * t82 / 0.2e1) * t105 + t116 * mrSges(5,3) + (-m(5) * t116 + t115) * t106) * qJD(4); (t16 * t44 + t17 * t43 + t89 * t92) * t167 - t39 * t42 - t112 * t20 - t38 * t41 - t113 * t19 + 0.2e1 * t89 * t40 + 0.2e1 * t92 * t18 - t104 * t73 + t101 * t72 + 0.2e1 * qJ(3) * t71 + (-t146 - t149) * qJD(4) + 0.2e1 * (mrSges(4,3) + t81 + (m(4) + m(5)) * qJ(3)) * qJD(3) + t107 * t135; t101 * t49 + t104 * t50 - t112 * t13 + t113 * t14 + t38 * t47 - t39 * t48 + t115 * qJD(4) + t133 * t141 + m(6) * t108 + m(5) * (-t116 * qJD(4) + t171); -m(6) * t107 + t124 * t135; -0.2e1 * m(6) * t124; -Ifges(5,5) * t131 + t12 * mrSges(5,1) - t11 * mrSges(5,2) + (m(6) * (t10 * t136 + t100 * t2 + t103 * t3 - t9 * t137) - t48 * t137 + t103 * t13 + t47 * t136 + t100 * t14) * pkin(4) + t111 + t172; ((-mrSges(5,2) * t106 - Ifges(5,6)) * t104 + (-mrSges(5,1) * t106 - Ifges(5,5)) * t101) * qJD(4) + (m(6) * (t100 * t16 + t103 * t17 + (-t100 * t43 + t103 * t44) * qJD(5)) + t169 * mrSges(6,3)) * pkin(4) + t114; -m(6) * t169 * pkin(4) - t81 * qJD(4) + t128; 0.2e1 * t63; t134 + t172; t114; t128; t63; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
