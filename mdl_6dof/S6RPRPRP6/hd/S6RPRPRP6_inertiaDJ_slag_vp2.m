% Calculate time derivative of joint inertia matrix for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:18:04
% EndTime: 2019-03-09 03:18:07
% DurationCPUTime: 2.07s
% Computational Cost: add. (2350->294), mult. (5177->407), div. (0->0), fcn. (4728->6), ass. (0->132)
t172 = Ifges(6,6) + Ifges(7,6);
t155 = cos(qJ(3));
t95 = sin(pkin(9));
t96 = cos(pkin(9));
t98 = sin(qJ(3));
t70 = t155 * t95 + t98 * t96;
t174 = t172 * t70;
t173 = Ifges(6,5) + Ifges(7,5);
t171 = -Ifges(6,3) - Ifges(7,3);
t99 = cos(qJ(5));
t151 = Ifges(7,4) * t99;
t97 = sin(qJ(5));
t110 = Ifges(7,1) * t97 + t151;
t153 = Ifges(6,4) * t99;
t111 = Ifges(6,1) * t97 + t153;
t120 = t155 * t96;
t148 = t95 * t98;
t69 = -t120 + t148;
t136 = t173 * t70 + (t110 + t111) * t69;
t170 = t136 * t97;
t152 = Ifges(7,4) * t97;
t108 = Ifges(7,2) * t99 + t152;
t154 = Ifges(6,4) * t97;
t109 = Ifges(6,2) * t99 + t154;
t137 = t174 + (t108 + t109) * t69;
t169 = t137 * t99;
t130 = qJD(5) * t99;
t66 = t70 * qJD(3);
t103 = t69 * t130 + t97 * t66;
t131 = qJD(5) * t97;
t126 = t69 * t131;
t146 = t99 * t66;
t102 = t126 - t146;
t168 = -t172 * t99 - t173 * t97;
t164 = -2 * mrSges(7,3);
t167 = 2 * mrSges(5,1) + 2 * mrSges(4,3);
t121 = -pkin(2) * t96 - pkin(1);
t104 = -qJ(4) * t70 + t121;
t158 = pkin(3) + pkin(8);
t29 = t158 * t69 + t104;
t140 = pkin(7) + qJ(2);
t79 = t140 * t95;
t80 = t140 * t96;
t48 = t155 * t79 + t80 * t98;
t39 = pkin(4) * t70 + t48;
t9 = t99 * t29 + t97 * t39;
t119 = qJD(3) * t155;
t65 = qJD(3) * t148 - t119 * t96;
t106 = qJ(4) * t65 - qJD(4) * t70;
t17 = t158 * t66 + t106;
t49 = t155 * t80 - t98 * t79;
t36 = qJD(2) * t70 + qJD(3) * t49;
t21 = -t65 * pkin(4) + t36;
t3 = t39 * t130 - t131 * t29 + t99 * t17 + t97 * t21;
t19 = t99 * t21;
t4 = -t9 * qJD(5) - t17 * t97 + t19;
t166 = t3 * t97 + t4 * t99;
t35 = (qJD(2) * t95 + qJD(3) * t80) * t98 - qJD(2) * t120 + t79 * t119;
t165 = 2 * m(7);
t28 = pkin(3) * t66 + t106;
t163 = -0.2e1 * t28;
t162 = m(7) * pkin(5);
t150 = t69 * t97;
t149 = t69 * t99;
t144 = mrSges(5,2) - mrSges(4,1);
t143 = mrSges(6,2) + mrSges(7,2);
t23 = mrSges(7,2) * t65 - mrSges(7,3) * t102;
t24 = mrSges(6,2) * t65 - mrSges(6,3) * t102;
t139 = t23 + t24;
t25 = -mrSges(7,1) * t65 - mrSges(7,3) * t103;
t26 = -mrSges(6,1) * t65 - mrSges(6,3) * t103;
t138 = t25 + t26;
t44 = mrSges(7,1) * t70 - mrSges(7,3) * t150;
t45 = mrSges(6,1) * t70 - mrSges(6,3) * t150;
t135 = -t44 - t45;
t46 = -mrSges(7,2) * t70 + mrSges(7,3) * t149;
t47 = -mrSges(6,2) * t70 + mrSges(6,3) * t149;
t134 = t46 + t47;
t132 = qJ(6) * t69;
t129 = t99 * qJD(6);
t128 = qJ(6) + t158;
t127 = 0.2e1 * t66;
t83 = -Ifges(7,2) * t97 + t151;
t84 = -Ifges(6,2) * t97 + t153;
t124 = t83 / 0.2e1 + t84 / 0.2e1;
t85 = Ifges(7,1) * t99 - t152;
t86 = Ifges(6,1) * t99 - t154;
t123 = t85 / 0.2e1 + t86 / 0.2e1;
t122 = mrSges(7,1) + t162;
t118 = -t29 - t132;
t78 = t128 * t99;
t92 = mrSges(7,1) * t130;
t71 = -mrSges(7,2) * t131 + t92;
t38 = t99 * t39;
t5 = pkin(5) * t70 + t118 * t97 + t38;
t6 = t132 * t99 + t9;
t114 = t5 * t97 - t6 * t99;
t8 = -t29 * t97 + t38;
t113 = t8 * t97 - t9 * t99;
t112 = mrSges(6,1) * t99 - mrSges(6,2) * t97;
t107 = -t35 * t49 + t36 * t48;
t14 = t102 * mrSges(7,1) + t103 * mrSges(7,2);
t101 = t173 * t103 + t172 * t146 + t171 * t65;
t20 = -pkin(4) * t66 - t35;
t91 = pkin(5) * t97 + qJ(4);
t87 = pkin(5) * t130 + qJD(4);
t82 = mrSges(6,1) * t97 + mrSges(6,2) * t99;
t81 = mrSges(7,1) * t97 + mrSges(7,2) * t99;
t77 = t128 * t97;
t76 = t111 * qJD(5);
t75 = t110 * qJD(5);
t74 = t109 * qJD(5);
t73 = t108 * qJD(5);
t72 = t112 * qJD(5);
t62 = t65 * mrSges(4,2);
t61 = t65 * mrSges(5,3);
t60 = -qJD(5) * t78 - t97 * qJD(6);
t59 = t128 * t131 - t129;
t43 = t112 * t69;
t42 = (-mrSges(7,1) * t99 + mrSges(7,2) * t97) * t69;
t41 = pkin(3) * t69 + t104;
t40 = -t69 * pkin(4) + t49;
t22 = (-pkin(5) * t99 - pkin(4)) * t69 + t49;
t15 = mrSges(6,1) * t102 + mrSges(6,2) * t103;
t13 = Ifges(6,1) * t103 - Ifges(6,4) * t102 - t65 * Ifges(6,5);
t12 = Ifges(7,1) * t103 - Ifges(7,4) * t102 - t65 * Ifges(7,5);
t11 = Ifges(6,4) * t103 - Ifges(6,2) * t102 - t65 * Ifges(6,6);
t10 = Ifges(7,4) * t103 - Ifges(7,2) * t102 - t65 * Ifges(7,6);
t7 = pkin(5) * t102 + t20;
t2 = -qJ(6) * t102 + t129 * t69 + t3;
t1 = -pkin(5) * t65 + t19 + t118 * t130 + (-qJ(6) * t66 - qJD(5) * t39 - qJD(6) * t69 - t17) * t97;
t16 = [0.2e1 * t40 * t15 + 0.2e1 * t7 * t42 - 0.2e1 * t20 * t43 + 0.2e1 * t1 * t44 + 0.2e1 * t4 * t45 + 0.2e1 * t2 * t46 + 0.2e1 * t3 * t47 + 0.2e1 * t6 * t23 + 0.2e1 * t9 * t24 + 0.2e1 * t5 * t25 + 0.2e1 * t8 * t26 + 0.2e1 * t22 * t14 + 0.2e1 * t41 * t61 - 0.2e1 * t121 * t62 - t65 * t48 * t167 + (0.2e1 * mrSges(4,1) * t121 - 0.2e1 * t41 * mrSges(5,2) - t167 * t49 + t169 + t170) * t66 + 0.2e1 * m(4) * t107 + 0.2e1 * m(5) * (t28 * t41 + t107) + 0.2e1 * m(6) * (t20 * t40 + t3 * t9 + t4 * t8) + (t1 * t5 + t2 * t6 + t22 * t7) * t165 + (mrSges(5,3) * t163 + (-Ifges(4,4) - Ifges(5,6)) * t127 + t36 * t167 + (-(2 * Ifges(4,1)) - (2 * Ifges(5,2)) + t171) * t65 + t101) * t70 + (mrSges(5,2) * t163 + (t10 + t11) * t99 + (t12 + t13) * t97 + (Ifges(5,3) + Ifges(4,2)) * t127 + t35 * t167 + (0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + t168) * t65 + (t136 * t99 + (-t137 - t174) * t97) * qJD(5)) * t69 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2) * (t95 ^ 2 + t96 ^ 2); t61 - t62 + t139 * t99 - t138 * t97 - t144 * t66 + (-t134 * t97 + t135 * t99) * qJD(5) + m(7) * (-t1 * t97 + t2 * t99 + (-t5 * t99 - t6 * t97) * qJD(5)) + m(6) * (t3 * t99 - t4 * t97 + (-t8 * t99 - t9 * t97) * qJD(5)) + m(5) * t28; 0; qJ(4) * t15 - qJD(4) * t43 + t91 * t14 + t20 * t82 + t22 * t71 - t77 * t23 - t78 * t25 + t40 * t72 + t87 * t42 + t59 * t44 + t60 * t46 + t7 * t81 + t144 * t36 + (-mrSges(5,3) + mrSges(4,2)) * t35 + (-t1 * mrSges(7,3) - t4 * mrSges(6,3) - t158 * t26 + t12 / 0.2e1 + t13 / 0.2e1) * t99 + (-t2 * mrSges(7,3) - t3 * mrSges(6,3) - t158 * t24 - t10 / 0.2e1 - t11 / 0.2e1) * t97 + (pkin(3) * mrSges(5,1) + Ifges(5,4) - Ifges(4,5) + (-Ifges(7,5) / 0.2e1 - Ifges(6,5) / 0.2e1) * t99 + (Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t97) * t65 + m(7) * (-t1 * t78 - t2 * t77 + t22 * t87 + t5 * t59 + t6 * t60 + t7 * t91) + m(5) * (-pkin(3) * t36 - qJ(4) * t35 + qJD(4) * t49) + m(6) * (qJ(4) * t20 + qJD(4) * t40 - t158 * t166) + (-qJD(4) * mrSges(5,1) + (-t73 / 0.2e1 - t74 / 0.2e1) * t99 + (-t75 / 0.2e1 - t76 / 0.2e1) * t97) * t69 + (-qJ(4) * mrSges(5,1) + t123 * t97 + t124 * t99 + Ifges(5,5) - Ifges(4,6)) * t66 + (t114 * mrSges(7,3) + t113 * mrSges(6,3) + (t123 * t99 - t124 * t97) * t69 - (-m(6) * t113 - t97 * t45 + t99 * t47) * t158 + t168 * t70 / 0.2e1 - t170 / 0.2e1 - t169 / 0.2e1) * qJD(5); m(7) * (-t59 * t97 + t60 * t99 + (t77 * t97 + t78 * t99) * qJD(5)); 0.2e1 * qJ(4) * t72 + 0.2e1 * t87 * t81 + 0.2e1 * t91 * t71 + (-t59 * t78 - t60 * t77 + t87 * t91) * t165 + (t59 * t164 - t75 - t76) * t99 + (t60 * t164 + t73 + t74) * t97 + 0.2e1 * (mrSges(5,3) + t82 + (m(5) + m(6)) * qJ(4)) * qJD(4) + ((-t164 * t77 - t83 - t84) * t99 + (t78 * t164 - t85 - t86) * t97) * qJD(5); -t65 * mrSges(5,1) + t138 * t99 + t139 * t97 + (t134 * t99 + t135 * t97) * qJD(5) + m(7) * (-qJD(5) * t114 + t1 * t99 + t2 * t97) + m(6) * (-qJD(5) * t113 + t166) + m(5) * t36; 0; m(7) * (t59 * t99 + t60 * t97 + (-t77 * t99 + t78 * t97) * qJD(5)); 0; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 - t172 * t126 + (m(7) * t1 + t25) * pkin(5) + t101; -t92 + ((-mrSges(6,1) - t162) * t99 + t143 * t97) * qJD(5); -t60 * mrSges(7,2) + t122 * t59 + ((mrSges(6,2) * t158 - t172) * t99 + (mrSges(6,1) * t158 + (mrSges(7,3) * pkin(5)) - t173) * t97) * qJD(5); (-t143 * t99 + (-mrSges(6,1) - t122) * t97) * qJD(5); 0; m(7) * t7 + t14; 0; m(7) * t87 + t71; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
