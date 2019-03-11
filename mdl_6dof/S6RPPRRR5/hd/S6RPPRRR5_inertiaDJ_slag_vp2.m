% Calculate time derivative of joint inertia matrix for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:16
% EndTime: 2019-03-09 02:28:20
% DurationCPUTime: 1.47s
% Computational Cost: add. (2127->241), mult. (4196->358), div. (0->0), fcn. (3522->6), ass. (0->127)
t82 = sin(qJ(6));
t84 = cos(qJ(6));
t60 = -mrSges(7,1) * t84 + mrSges(7,2) * t82;
t186 = mrSges(6,1) - t60;
t76 = t82 ^ 2;
t78 = t84 ^ 2;
t184 = t76 + t78;
t85 = cos(qJ(4));
t142 = qJD(4) * t85;
t167 = sin(qJ(5));
t116 = qJD(5) * t167;
t182 = -qJD(4) * t167 - t116;
t168 = cos(qJ(5));
t114 = t168 * qJD(5);
t183 = t168 * qJD(4) + t114;
t83 = sin(qJ(4));
t39 = t182 * t83 + t183 * t85;
t141 = qJD(6) * t82;
t125 = t168 * t85;
t54 = t167 * t83 - t125;
t132 = t54 * t141;
t38 = t182 * t85 - t183 * t83;
t157 = t38 * t84;
t93 = t132 + t157;
t10 = mrSges(7,1) * t39 - mrSges(7,3) * t93;
t140 = qJD(6) * t84;
t158 = t38 * t82;
t94 = t140 * t54 - t158;
t11 = -mrSges(7,2) * t39 + mrSges(7,3) * t94;
t104 = -t82 * t10 + t84 * t11;
t153 = t82 * mrSges(7,3);
t121 = t167 * t85;
t55 = t168 * t83 + t121;
t32 = -mrSges(7,2) * t55 + t153 * t54;
t155 = t54 * t84;
t33 = mrSges(7,1) * t55 + mrSges(7,3) * t155;
t185 = -t33 * t140 - t32 * t141 + t104;
t181 = -t167 * t39 - t168 * t38;
t81 = pkin(1) + qJ(3);
t70 = t83 * pkin(4) + t81;
t27 = pkin(5) * t55 + pkin(9) * t54 + t70;
t80 = qJ(2) - pkin(7);
t169 = pkin(8) - t80;
t59 = t169 * t83;
t43 = -t121 * t169 - t168 * t59;
t15 = t27 * t84 - t43 * t82;
t16 = t27 * t82 + t43 * t84;
t179 = -t15 * t140 - t16 * t141;
t143 = qJD(4) * t83;
t178 = -mrSges(5,1) * t143 - mrSges(5,2) * t142;
t102 = t84 * t32 - t82 * t33;
t77 = t83 ^ 2;
t79 = t85 ^ 2;
t177 = -2 * mrSges(6,3);
t175 = m(6) * pkin(4);
t108 = mrSges(7,1) * t82 + mrSges(7,2) * t84;
t56 = t108 * qJD(6);
t172 = pkin(5) * t56;
t42 = t125 * t169 - t167 * t59;
t139 = t83 * qJD(2);
t49 = -t142 * t169 + t139;
t138 = t85 * qJD(2);
t89 = t143 * t169 + t138;
t13 = -qJD(5) * t42 + t167 * t89 + t168 * t49;
t63 = pkin(4) * t142 + qJD(3);
t17 = t39 * pkin(5) - t38 * pkin(9) + t63;
t2 = qJD(6) * t15 + t13 * t84 + t17 * t82;
t171 = t2 * t84;
t3 = -qJD(6) * t16 - t13 * t82 + t17 * t84;
t170 = t3 * t82;
t165 = Ifges(7,4) * t82;
t164 = Ifges(7,4) * t84;
t163 = Ifges(7,6) * t82;
t14 = qJD(5) * t43 + t167 * t49 - t168 * t89;
t162 = t14 * t42;
t161 = t14 * t54;
t160 = t38 * t42;
t159 = t38 * t54;
t156 = t54 * mrSges(6,3);
t73 = -pkin(4) * t168 - pkin(5);
t154 = t73 * t56;
t147 = Ifges(7,5) * t157 + Ifges(7,3) * t39;
t146 = t77 + t79;
t145 = pkin(4) * qJD(5);
t144 = qJD(3) * t81;
t137 = m(7) * t145;
t128 = t82 * t168;
t127 = t84 * t168;
t123 = t167 * t42;
t122 = t167 * t54;
t120 = t184 * t39;
t118 = -(2 * Ifges(6,4)) - t163;
t113 = t146 * qJD(2);
t112 = pkin(4) * t114;
t111 = pkin(4) * t116;
t110 = t85 * mrSges(5,1) - t83 * mrSges(5,2);
t109 = t39 * mrSges(6,1) + t38 * mrSges(6,2);
t107 = Ifges(7,1) * t84 - t165;
t106 = -Ifges(7,2) * t82 + t164;
t105 = Ifges(7,5) * t82 + Ifges(7,6) * t84;
t103 = t15 * t82 - t16 * t84;
t101 = -t39 * t43 + t160;
t100 = mrSges(7,3) * t112;
t99 = mrSges(6,2) * t112;
t98 = mrSges(6,1) * t111;
t97 = t60 * t111;
t96 = t76 * t100;
t95 = t78 * t100;
t57 = t106 * qJD(6);
t58 = t107 * qJD(6);
t61 = Ifges(7,2) * t84 + t165;
t62 = Ifges(7,1) * t82 + t164;
t92 = t62 * t140 - t141 * t61 + t84 * t57 + t82 * t58;
t91 = t186 * t38 + t54 * t56 + (mrSges(7,3) * t184 - mrSges(6,2)) * t39;
t90 = t184 * t168;
t88 = -t170 + (-t15 * t84 - t16 * t82) * qJD(6);
t87 = t88 + t171;
t20 = Ifges(7,6) * t55 - t106 * t54;
t21 = Ifges(7,5) * t55 - t107 * t54;
t6 = Ifges(7,4) * t93 + Ifges(7,2) * t94 + Ifges(7,6) * t39;
t7 = Ifges(7,1) * t93 + Ifges(7,4) * t94 + Ifges(7,5) * t39;
t74 = Ifges(7,5) * t140;
t86 = -t13 * mrSges(6,2) + mrSges(7,3) * t171 - t20 * t141 / 0.2e1 - t61 * t158 / 0.2e1 + t62 * t157 / 0.2e1 + t42 * t56 + Ifges(6,5) * t38 - t58 * t155 / 0.2e1 + t55 * (-Ifges(7,6) * t141 + t74) / 0.2e1 + t84 * t6 / 0.2e1 + (t105 / 0.2e1 - Ifges(6,6)) * t39 - t186 * t14 + (t54 * t61 + t21) * t140 / 0.2e1 + (t7 + (qJD(6) * t62 + t57) * t54) * t82 / 0.2e1;
t72 = pkin(4) * t167 + pkin(9);
t28 = t108 * t54;
t8 = -mrSges(7,1) * t94 + mrSges(7,2) * t93;
t1 = [t21 * t157 - t20 * t158 + 0.2e1 * t70 * t109 + 0.2e1 * t42 * t8 - 0.2e1 * t14 * t28 + 0.2e1 * t2 * t32 + 0.2e1 * t3 * t33 + 0.2e1 * t15 * t10 + 0.2e1 * t16 * t11 + 0.2e1 * t101 * mrSges(6,3) + 0.2e1 * m(5) * (t113 * t80 + t144) + 0.2e1 * m(4) * (qJ(2) * qJD(2) + t144) + 0.2e1 * m(6) * (t13 * t43 + t63 * t70 + t162) + 0.2e1 * m(7) * (t15 * t3 + t16 * t2 + t162) + (0.2e1 * t63 * mrSges(6,1) + t13 * t177 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t39 + t118 * t38 + t147) * t55 + (-0.2e1 * t63 * mrSges(6,2) + t14 * t177 - 0.2e1 * Ifges(6,1) * t38 + t82 * t6 - t84 * t7 + (-Ifges(7,5) * t84 - t118) * t39 + (t105 * t55 + t84 * t20 + t82 * t21) * qJD(6)) * t54 + 0.2e1 * (t81 * t110 + (-t79 + t77) * Ifges(5,4)) * qJD(4) + 0.2e1 * (m(3) * qJ(2) - mrSges(5,3) * t146 + mrSges(4,2) + mrSges(3,3)) * qJD(2) + 0.2e1 * (-Ifges(5,1) + Ifges(5,2)) * t83 * t142 + 0.2e1 * (t83 * mrSges(5,1) + t85 * mrSges(5,2) + mrSges(4,3)) * qJD(3); -t84 * t10 - t82 * t11 - t102 * qJD(6) - t110 * qJD(4) + (-m(5) - m(4)) * qJD(3) + m(7) * (qJD(6) * t103 - t2 * t82 - t3 * t84) - m(6) * t63 - t109; 0; m(4) * qJD(2) + t54 * t8 + t102 * t39 + (t28 + 0.2e1 * t156) * t38 + (t39 * t177 + (-t32 * t82 - t33 * t84) * qJD(6) + t104) * t55 + m(7) * (-t103 * t39 + t55 * t87 - t160 + t161) + m(6) * (t13 * t55 - t101 + t161) + m(5) * t113; 0; 0.2e1 * m(6) * (t39 * t55 - t159) + 0.2e1 * m(7) * (t120 * t55 - t159); mrSges(5,1) * t138 - Ifges(5,6) * t142 - mrSges(5,2) * t139 - Ifges(5,5) * t143 + t73 * t8 + t86 - t3 * t153 + (-t168 * t14 + t167 * t13 + (t168 * t43 + t123) * qJD(5)) * t175 + m(7) * (t73 * t14 + (t127 * t16 - t128 * t15 + t123) * t145) + t178 * t80 + t181 * mrSges(6,3) * pkin(4) + (-t28 - t156) * t111 + t179 * mrSges(7,3) + (-mrSges(6,3) * t55 + t102) * t112 + (m(7) * t87 + t185) * t72; 0; ((t168 * t55 + t122) * qJD(5) - t181) * t175 + m(7) * (-t73 * t38 + t72 * t120 + (t55 * t90 + t122) * t145) + t91 + t178; 0.2e1 * t95 + 0.2e1 * t96 - 0.2e1 * t99 - 0.2e1 * t98 + 0.2e1 * t97 + 0.2e1 * t154 + 0.2e1 * (t167 * t73 + t72 * t90) * t137 + t92; t86 + t88 * mrSges(7,3) + (-m(7) * t14 - t8) * pkin(5) + (m(7) * (-t170 + t171 + t179) + t185) * pkin(9); 0; m(7) * (pkin(5) * t38 + pkin(9) * t120) + t91; t95 + t96 - t99 - t98 + t97 + t154 + (-pkin(5) * t167 + pkin(9) * t90) * t137 - t172 + t92; t92 - 0.2e1 * t172; mrSges(7,1) * t3 - mrSges(7,2) * t2 + Ifges(7,5) * t132 + Ifges(7,6) * t94 + t147; t56; (t141 * t55 - t39 * t84) * mrSges(7,2) + (-t140 * t55 - t39 * t82) * mrSges(7,1); t74 + (-mrSges(7,1) * t128 - mrSges(7,2) * t127) * t145 + (t60 * t72 - t163) * qJD(6); t74 + (pkin(9) * t60 - t163) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
