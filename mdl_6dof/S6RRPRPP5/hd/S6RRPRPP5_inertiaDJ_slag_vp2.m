% Calculate time derivative of joint inertia matrix for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:32
% EndTime: 2019-03-09 10:03:38
% DurationCPUTime: 2.45s
% Computational Cost: add. (1569->387), mult. (3487->515), div. (0->0), fcn. (2249->4), ass. (0->166)
t183 = pkin(3) + pkin(7);
t193 = Ifges(6,4) + Ifges(5,5);
t192 = Ifges(6,2) + Ifges(5,3);
t110 = sin(qJ(4));
t112 = cos(qJ(4));
t113 = cos(qJ(2));
t153 = qJD(4) * t113;
t111 = sin(qJ(2));
t157 = qJD(2) * t111;
t117 = t110 * t153 + t112 * t157;
t143 = t110 * t157;
t144 = t112 * t153;
t116 = t143 - t144;
t191 = 2 * qJ(3);
t164 = qJ(5) * t112;
t182 = pkin(4) + pkin(5);
t190 = t110 * t182 - t164;
t165 = qJ(5) * t110;
t119 = -t112 * t182 - t165;
t189 = 2 * m(7);
t188 = -2 * pkin(1);
t187 = 2 * mrSges(7,3);
t139 = pkin(2) * t157 - qJD(3) * t111;
t156 = qJD(2) * t113;
t48 = -qJ(3) * t156 + t139;
t186 = 0.2e1 * t48;
t141 = -qJ(3) * t111 - pkin(1);
t78 = -pkin(2) * t113 + t141;
t185 = -0.2e1 * t78;
t184 = m(6) + m(7);
t115 = -pkin(2) - pkin(8);
t181 = -mrSges(6,1) - mrSges(7,1);
t180 = -mrSges(6,2) + mrSges(7,3);
t179 = mrSges(7,2) + mrSges(6,3);
t178 = Ifges(3,4) + Ifges(4,6);
t177 = -Ifges(6,6) + Ifges(7,6);
t29 = -mrSges(5,2) * t156 + mrSges(5,3) * t117;
t33 = mrSges(6,2) * t117 + mrSges(6,3) * t156;
t176 = t29 + t33;
t31 = mrSges(5,1) * t156 - mrSges(5,3) * t116;
t154 = qJD(4) * t112;
t91 = mrSges(6,2) * t143;
t32 = t91 + (-mrSges(6,1) * qJD(2) - mrSges(6,2) * t154) * t113;
t175 = t31 - t32;
t54 = t113 * t115 + t141;
t88 = t183 * t111;
t23 = t110 * t88 + t112 * t54;
t163 = t110 * t113;
t68 = mrSges(5,1) * t111 + mrSges(5,3) * t163;
t69 = -mrSges(6,1) * t111 - mrSges(6,2) * t163;
t174 = -t68 + t69;
t161 = t112 * t113;
t71 = -mrSges(5,2) * t111 - mrSges(5,3) * t161;
t72 = -mrSges(6,2) * t161 + mrSges(6,3) * t111;
t173 = t71 + t72;
t172 = Ifges(5,4) * t110;
t171 = Ifges(5,4) * t112;
t170 = Ifges(7,4) * t110;
t169 = Ifges(7,4) * t112;
t168 = Ifges(6,5) * t110;
t167 = Ifges(6,5) * t112;
t166 = t111 * Ifges(7,5);
t162 = t110 * t115;
t160 = t112 * t115;
t159 = qJ(6) + t115;
t106 = qJD(5) * t110;
t158 = qJ(5) * t154 + t106;
t89 = t183 * t113;
t155 = qJD(4) * t110;
t152 = t110 * qJD(6);
t151 = t112 * qJD(6);
t150 = Ifges(7,5) - t193;
t129 = Ifges(7,1) * t110 - t169;
t40 = -t113 * t129 - t166;
t130 = Ifges(6,1) * t110 - t167;
t41 = t111 * Ifges(6,4) - t113 * t130;
t131 = Ifges(5,1) * t110 + t171;
t42 = t111 * Ifges(5,5) - t113 * t131;
t149 = -t40 - t41 - t42;
t17 = t111 * qJ(5) + t23;
t148 = t110 * t111 * mrSges(7,3);
t147 = -Ifges(5,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t146 = m(4) * pkin(7) + mrSges(4,1);
t22 = -t110 * t54 + t112 * t88;
t140 = qJD(4) * t159;
t138 = -qJD(5) * t112 + qJD(3);
t137 = -Ifges(5,5) / 0.2e1 + Ifges(7,5) / 0.2e1 - Ifges(6,4) / 0.2e1;
t82 = Ifges(6,3) * t110 + t167;
t83 = Ifges(7,2) * t110 + t169;
t84 = -Ifges(5,2) * t110 + t171;
t136 = -t82 / 0.2e1 - t83 / 0.2e1 + t84 / 0.2e1;
t85 = Ifges(7,1) * t112 + t170;
t86 = Ifges(6,1) * t112 + t168;
t87 = Ifges(5,1) * t112 - t172;
t135 = t86 / 0.2e1 + t87 / 0.2e1 + t85 / 0.2e1;
t134 = mrSges(5,1) * t112 - mrSges(5,2) * t110;
t81 = t110 * mrSges(5,1) + t112 * mrSges(5,2);
t133 = mrSges(6,1) * t112 + mrSges(6,3) * t110;
t79 = t110 * mrSges(6,1) - t112 * mrSges(6,3);
t132 = -mrSges(7,1) * t112 - mrSges(7,2) * t110;
t128 = Ifges(5,2) * t112 + t172;
t127 = -Ifges(7,2) * t112 + t170;
t126 = -Ifges(6,3) * t112 + t168;
t125 = pkin(4) * t112 + t165;
t124 = -pkin(4) * t110 + t164;
t18 = -pkin(4) * t111 - t22;
t123 = t110 * t18 + t112 * t17;
t122 = -t110 * t22 + t112 * t23;
t35 = (pkin(8) * t111 - qJ(3) * t113) * qJD(2) + t139;
t74 = t183 * t156;
t6 = -t110 * t35 + t112 * t74 - t54 * t154 - t88 * t155;
t121 = Ifges(5,6) * t117 + t143 * t193 + t156 * t192;
t120 = t110 * (m(6) * t115 + t180);
t5 = t110 * t74 + t112 * t35 + t88 * t154 - t155 * t54;
t37 = Ifges(6,6) * t111 - t113 * t126;
t38 = -t111 * Ifges(7,6) - t113 * t127;
t39 = t111 * Ifges(5,6) - t113 * t128;
t118 = t111 * t177 - t37 - t38 + t39;
t59 = t132 * qJD(4);
t20 = mrSges(7,1) * t117 + mrSges(7,2) * t116;
t3 = qJ(5) * t156 + t111 * qJD(5) + t5;
t102 = Ifges(6,6) * t154;
t94 = mrSges(7,3) * t144;
t80 = -mrSges(7,1) * t110 + mrSges(7,2) * t112;
t77 = qJ(3) - t124;
t76 = t159 * t112;
t75 = t159 * t110;
t73 = t183 * t157;
t70 = mrSges(7,2) * t111 + mrSges(7,3) * t161;
t67 = -mrSges(7,1) * t111 + mrSges(7,3) * t163;
t66 = t131 * qJD(4);
t65 = t130 * qJD(4);
t64 = t129 * qJD(4);
t63 = t128 * qJD(4);
t62 = t127 * qJD(4);
t61 = t126 * qJD(4);
t60 = t134 * qJD(4);
t58 = t133 * qJD(4);
t53 = t134 * t113;
t52 = t132 * t113;
t51 = t133 * t113;
t50 = -qJ(3) - t190;
t45 = qJD(4) * t125 + t138;
t44 = t112 * t140 + t152;
t43 = t110 * t140 - t151;
t34 = t113 * t125 + t89;
t30 = t94 + (-mrSges(7,1) * t113 - t148) * qJD(2);
t28 = mrSges(7,2) * t156 - mrSges(7,3) * t117;
t27 = qJD(4) * t119 - t138;
t24 = t113 * t119 - t89;
t21 = -mrSges(5,1) * t117 + mrSges(5,2) * t116;
t19 = -mrSges(6,1) * t117 - mrSges(6,3) * t116;
t16 = -t87 * t153 + (t113 * Ifges(5,5) + t111 * t131) * qJD(2);
t15 = -t86 * t153 + (t113 * Ifges(6,4) + t111 * t130) * qJD(2);
t14 = -t85 * t153 + (-t113 * Ifges(7,5) + t111 * t129) * qJD(2);
t13 = -t84 * t153 + (t113 * Ifges(5,6) + t111 * t128) * qJD(2);
t12 = -t83 * t153 + (-t113 * Ifges(7,6) + t111 * t127) * qJD(2);
t11 = -t82 * t153 + (t113 * Ifges(6,6) + t111 * t126) * qJD(2);
t10 = qJ(6) * t161 + t17;
t9 = qJ(6) * t163 - t111 * t182 - t22;
t8 = (qJD(4) * t124 + t106) * t113 + (-t125 - t183) * t157;
t7 = (qJD(4) * t190 - t106) * t113 + (-t119 + t183) * t157;
t4 = -pkin(4) * t156 - t6;
t2 = -qJ(6) * t117 + t113 * t151 + t3;
t1 = -qJ(6) * t143 + (qJ(6) * t154 - qJD(2) * t182 + t152) * t113 - t6;
t25 = [0.2e1 * t89 * t21 + 0.2e1 * t1 * t67 + 0.2e1 * t6 * t68 + 0.2e1 * t4 * t69 + 0.2e1 * t2 * t70 + 0.2e1 * t5 * t71 + 0.2e1 * t3 * t72 - 0.2e1 * t73 * t53 + 0.2e1 * t8 * t51 + 0.2e1 * t7 * t52 + 0.2e1 * t10 * t28 + 0.2e1 * t23 * t29 + 0.2e1 * t9 * t30 + 0.2e1 * t22 * t31 + 0.2e1 * t18 * t32 + 0.2e1 * t17 * t33 + 0.2e1 * t34 * t19 + m(4) * t78 * t186 + 0.2e1 * t24 * t20 + 0.2e1 * m(5) * (t22 * t6 + t23 * t5 - t73 * t89) + (t1 * t9 + t10 * t2 + t24 * t7) * t189 + 0.2e1 * m(6) * (t17 * t3 + t18 * t4 + t34 * t8) + (-0.2e1 * t48 * mrSges(4,3) + ((mrSges(3,1) * t188) + mrSges(4,2) * t185 - 0.2e1 * t178 * t111 + (-t149 - t166) * t110 + t118 * t112) * qJD(2) + t121) * t111 + (mrSges(4,2) * t186 + (t11 + t12 - t13) * t112 + (-t14 - t15 - t16) * t110 + (t118 * t110 + (t111 * t150 + t149) * t112) * qJD(4) + (0.2e1 * t178 * t113 + mrSges(4,3) * t185 + (mrSges(3,2) * t188) + (-Ifges(5,6) - t177) * t161 + t150 * t163 + ((2 * Ifges(7,3)) + (2 * Ifges(4,2)) + (2 * Ifges(3,1)) - (2 * Ifges(4,3)) - (2 * Ifges(3,2)) + t192) * t111) * qJD(2)) * t113; t111 * t102 / 0.2e1 + t8 * t79 + t7 * t80 - t73 * t81 + t89 * t60 + t34 * t58 + t24 * t59 + t43 * t67 + t44 * t70 + t75 * t28 - t76 * t30 + t77 * t19 + t50 * t20 + t45 * t51 + t27 * t52 + qJD(3) * t53 + qJ(3) * t21 + m(7) * (-t1 * t76 + t10 * t44 + t2 * t75 + t24 * t27 + t43 * t9 + t50 * t7) + (t146 * qJD(3) + (-t61 / 0.2e1 - t62 / 0.2e1 + t63 / 0.2e1) * t112 + (t64 / 0.2e1 + t65 / 0.2e1 + t66 / 0.2e1) * t110) * t113 + m(5) * (-qJ(3) * t73 + qJD(3) * t89 + t160 * t6 + t162 * t5) + ((t37 / 0.2e1 + t38 / 0.2e1 - t39 / 0.2e1 - t23 * mrSges(5,3) + t10 * mrSges(7,3) - t17 * mrSges(6,2) + t147 * t111 - t135 * t113) * t112 + (-t40 / 0.2e1 - t41 / 0.2e1 - t42 / 0.2e1 + t9 * mrSges(7,3) - t18 * mrSges(6,2) + t22 * mrSges(5,3) + t136 * t113 + t137 * t111) * t110 + (m(5) * t122 + m(6) * t123 + t110 * t174 + t112 * t173) * t115) * qJD(4) + (t14 / 0.2e1 + t15 / 0.2e1 + t16 / 0.2e1 - t1 * mrSges(7,3) + t4 * mrSges(6,2) - t6 * mrSges(5,3) + t175 * t115) * t112 + (t11 / 0.2e1 + t12 / 0.2e1 - t13 / 0.2e1 + t2 * mrSges(7,3) - t3 * mrSges(6,2) - t5 * mrSges(5,3) + t176 * t115) * t110 + m(6) * (-t160 * t4 + t162 * t3 + t45 * t34 + t77 * t8) + ((-pkin(2) * mrSges(4,1) - Ifges(4,4) + Ifges(3,5) - t137 * t112 + (Ifges(6,6) / 0.2e1 + t147) * t110 + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2)) * pkin(7)) * t113 + (-qJ(3) * mrSges(4,1) + Ifges(4,5) - Ifges(3,6) + t136 * t112 + t135 * t110 + (-m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3)) * pkin(7)) * t111) * qJD(2); 0.2e1 * t27 * t80 + 0.2e1 * t50 * t59 + t60 * t191 + 0.2e1 * t77 * t58 + (t27 * t50 - t43 * t76 + t44 * t75) * t189 + (-0.2e1 * mrSges(7,3) * t43 - t64 - t65 - t66) * t112 + (t44 * t187 - t61 - t62 + t63) * t110 + ((t187 * t75 + t82 + t83 - t84) * t112 + (-t187 * t76 - t85 - t86 - t87) * t110) * qJD(4) + 0.2e1 * (m(6) * t77 + t79) * t45 + (0.2e1 * mrSges(4,3) + 0.2e1 * t81 + (m(4) + m(5)) * t191) * qJD(3); t146 * t156 + (-t30 + t175) * t112 + (t28 + t176) * t110 + ((t70 + t173) * t112 + (t67 + t174) * t110) * qJD(4) + m(7) * (-t1 * t112 + t110 * t2 + (t10 * t112 + t110 * t9) * qJD(4)) + m(6) * (qJD(4) * t123 + t110 * t3 - t112 * t4) + m(5) * (qJD(4) * t122 + t110 * t5 + t112 * t6); m(7) * (t110 * t44 - t112 * t43 + (-t110 * t76 + t112 * t75) * qJD(4)); 0; -t182 * t30 - pkin(4) * t32 + t121 - t1 * mrSges(7,1) + t2 * mrSges(7,2) + t3 * mrSges(6,3) - t4 * mrSges(6,1) - t5 * mrSges(5,2) + t6 * mrSges(5,1) + (Ifges(7,3) * qJD(2) + (t110 * t177 + t112 * t150) * qJD(4)) * t113 + (t70 + t72) * qJD(5) + (t28 + t33) * qJ(5) + m(7) * (qJ(5) * t2 + qJD(5) * t10 - t1 * t182) + m(6) * (-pkin(4) * t4 + qJ(5) * t3 + qJD(5) * t17) + (-Ifges(7,5) * t110 + t112 * t177) * t157; t102 + m(7) * (qJ(5) * t44 + qJD(5) * t75 - t182 * t43) - t43 * mrSges(7,1) + t44 * mrSges(7,2) + qJD(5) * t120 + ((qJ(5) * t180 - Ifges(5,6) - Ifges(7,6)) * t112 + (pkin(4) * mrSges(6,2) - mrSges(7,3) * t182 + t150) * t110 + (m(6) * t124 - t79 - t81) * t115) * qJD(4); (-mrSges(5,1) + t181) * t155 + m(6) * (-pkin(4) * t155 + t158) + m(7) * (-t155 * t182 + t158) + (-mrSges(5,2) + t179) * t154; 0.2e1 * (qJ(5) * t184 + t179) * qJD(5); -mrSges(6,2) * t144 + t91 + t94 + m(6) * t4 + m(7) * t1 + (t113 * t181 - t148) * qJD(2); m(7) * t43 + qJD(4) * t120; t184 * t155; 0; 0; m(7) * t7 + t20; m(7) * t27 + t59; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
