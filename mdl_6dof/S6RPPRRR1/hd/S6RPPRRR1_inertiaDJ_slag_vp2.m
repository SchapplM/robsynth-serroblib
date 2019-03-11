% Calculate time derivative of joint inertia matrix for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:17:58
% EndTime: 2019-03-09 02:18:02
% DurationCPUTime: 1.94s
% Computational Cost: add. (4677->236), mult. (9641->363), div. (0->0), fcn. (9859->10), ass. (0->113)
t100 = cos(qJ(6));
t97 = sin(qJ(6));
t143 = t100 ^ 2 + t97 ^ 2;
t101 = cos(qJ(5));
t102 = cos(qJ(4));
t95 = cos(pkin(11));
t134 = t102 * t95;
t94 = sin(pkin(11));
t99 = sin(qJ(4));
t77 = -t99 * t94 + t134;
t78 = t102 * t94 + t99 * t95;
t98 = sin(qJ(5));
t113 = t101 * t77 - t78 * t98;
t70 = t77 * qJD(4);
t71 = t78 * qJD(4);
t42 = qJD(5) * t113 + t101 * t70 - t71 * t98;
t177 = t143 * t42;
t171 = t143 * t101;
t86 = sin(pkin(10)) * pkin(1) + qJ(3);
t162 = pkin(7) + t86;
t74 = t162 * t94;
t75 = t162 * t95;
t59 = t102 * t75 - t99 * t74;
t48 = -t78 * qJD(3) - qJD(4) * t59;
t105 = -t70 * pkin(8) + t48;
t66 = t102 * t74;
t58 = -t75 * t99 - t66;
t110 = -pkin(8) * t78 + t58;
t50 = pkin(8) * t77 + t59;
t21 = -t101 * t110 + t98 * t50;
t47 = -qJD(4) * t66 + qJD(3) * t134 + (-qJD(3) * t94 - qJD(4) * t75) * t99;
t46 = -pkin(8) * t71 + t47;
t10 = -qJD(5) * t21 + t101 * t46 + t105 * t98;
t22 = t101 * t50 + t110 * t98;
t61 = t101 * t78 + t77 * t98;
t112 = -cos(pkin(10)) * pkin(1) - pkin(3) * t95 - pkin(2);
t62 = -pkin(4) * t77 + t112;
t25 = -pkin(5) * t113 - pkin(9) * t61 + t62;
t13 = t100 * t25 - t22 * t97;
t163 = t71 * pkin(4);
t43 = qJD(5) * t61 + t101 * t71 + t70 * t98;
t165 = pkin(5) * t43;
t17 = -pkin(9) * t42 + t163 + t165;
t2 = qJD(6) * t13 + t10 * t100 + t17 * t97;
t14 = t100 * t22 + t25 * t97;
t3 = -qJD(6) * t14 - t10 * t97 + t100 * t17;
t176 = t100 * t2 - t3 * t97;
t175 = (-t100 * t13 - t14 * t97) * qJD(6) + t176;
t132 = qJD(6) * t100;
t174 = t61 * t132 + t97 * t42;
t133 = qJD(6) * t97;
t130 = t61 * t133;
t139 = t100 * t42;
t109 = t130 - t139;
t116 = t100 * t14 - t13 * t97;
t170 = 2 * m(7);
t169 = -2 * mrSges(6,3);
t11 = qJD(5) * t22 - t101 * t105 + t98 * t46;
t168 = 0.2e1 * t11;
t167 = m(6) / 0.2e1;
t161 = Ifges(7,4) * t97;
t160 = Ifges(7,6) * t97;
t158 = t11 * t21;
t153 = t113 * t43;
t152 = t113 * t98;
t151 = t61 * t97;
t150 = t70 * t78;
t149 = t77 * t71;
t147 = t98 * mrSges(6,1);
t82 = -mrSges(7,1) * t100 + mrSges(7,2) * t97;
t146 = t98 * t82;
t145 = Ifges(7,5) * t139 + Ifges(7,3) * t43;
t144 = t43 * mrSges(6,1) + t42 * mrSges(6,2);
t142 = Ifges(7,4) * t100;
t141 = pkin(4) * qJD(5);
t138 = t100 * t61;
t137 = t101 * mrSges(6,2);
t131 = 0.2e1 * t163;
t12 = -mrSges(7,1) * t174 + mrSges(7,2) * t109;
t128 = m(7) * t11 - t12;
t68 = t70 * mrSges(5,2);
t125 = t71 * mrSges(5,1) + t68;
t124 = -(2 * Ifges(6,4)) - t160;
t122 = mrSges(7,3) * t171;
t121 = -t11 * t113 + t21 * t43;
t120 = mrSges(7,1) * t97 + mrSges(7,2) * t100;
t119 = Ifges(7,1) * t100 - t161;
t118 = -Ifges(7,2) * t97 + t142;
t117 = Ifges(7,5) * t97 + Ifges(7,6) * t100;
t15 = mrSges(7,1) * t43 + mrSges(7,3) * t109;
t16 = -mrSges(7,2) * t43 - mrSges(7,3) * t174;
t115 = t100 * t16 - t97 * t15;
t44 = mrSges(7,2) * t113 - mrSges(7,3) * t151;
t45 = -mrSges(7,1) * t113 - mrSges(7,3) * t138;
t114 = t100 * t44 - t97 * t45;
t79 = t120 * qJD(6);
t111 = mrSges(7,3) * t177 - t113 * t79 + t43 * t82 - t144;
t80 = t118 * qJD(6);
t81 = t119 * qJD(6);
t83 = Ifges(7,2) * t100 + t161;
t84 = Ifges(7,1) * t97 + t142;
t107 = t100 * t80 + t84 * t132 - t133 * t83 + t97 * t81;
t104 = -t45 * t132 - t44 * t133 + m(7) * (-t13 * t132 - t133 * t14 + t176) + t115;
t23 = -Ifges(7,6) * t113 + t118 * t61;
t24 = -Ifges(7,5) * t113 + t119 * t61;
t6 = -Ifges(7,4) * t109 - Ifges(7,2) * t174 + t43 * Ifges(7,6);
t7 = -Ifges(7,1) * t109 - Ifges(7,4) * t174 + t43 * Ifges(7,5);
t89 = Ifges(7,5) * t132;
t103 = -t10 * mrSges(6,2) + t21 * t79 + t24 * t132 / 0.2e1 + t84 * t139 / 0.2e1 + t97 * t7 / 0.2e1 + Ifges(6,5) * t42 + t100 * t6 / 0.2e1 - t80 * t151 / 0.2e1 + t81 * t138 / 0.2e1 - t113 * (-Ifges(7,6) * t133 + t89) / 0.2e1 + (t117 / 0.2e1 - Ifges(6,6)) * t43 - t174 * t83 / 0.2e1 - (t61 * t84 + t23) * t133 / 0.2e1 + (t82 - mrSges(6,1)) * t11 + t175 * mrSges(7,3);
t88 = -pkin(4) * t101 - pkin(5);
t87 = pkin(4) * t98 + pkin(9);
t35 = t120 * t61;
t1 = [0.2e1 * t112 * t125 + 0.2e1 * Ifges(5,1) * t150 - 0.2e1 * Ifges(5,2) * t149 + 0.2e1 * t62 * t144 + 0.2e1 * t2 * t44 + 0.2e1 * t3 * t45 + t35 * t168 - 0.2e1 * t21 * t12 + 0.2e1 * t13 * t15 + 0.2e1 * t14 * t16 + t22 * t43 * t169 + (0.2e1 * t21 * mrSges(6,3) + t100 * t24 - t97 * t23) * t42 + (t13 * t3 + t14 * t2 + t158) * t170 + 0.2e1 * m(5) * (t47 * t59 + t48 * t58) + 0.2e1 * m(6) * (t10 * t22 + t163 * t62 + t158) - (mrSges(6,1) * t131 + t10 * t169 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t43 + t124 * t42 + t145) * t113 + (mrSges(6,2) * t131 + mrSges(6,3) * t168 + 0.2e1 * Ifges(6,1) * t42 + t100 * t7 - t97 * t6 + (Ifges(7,5) * t100 + t124) * t43 + (-t100 * t23 + t113 * t117 - t97 * t24) * qJD(6)) * t61 + 0.2e1 * (t70 * t77 - t78 * t71) * Ifges(5,4) + 0.2e1 * (t47 * t77 - t48 * t78 - t58 * t70 - t59 * t71) * mrSges(5,3) + 0.2e1 * (m(4) * t86 + mrSges(4,3)) * qJD(3) * (t94 ^ 2 + t95 ^ 2); t113 * t12 + t43 * t35 + t114 * t42 + ((-t100 * t45 - t44 * t97) * qJD(6) + t115) * t61 + m(7) * (t116 * t42 + t175 * t61 + t121) + m(6) * (t10 * t61 + t22 * t42 + t121) + m(5) * (t47 * t78 + t48 * t77 - t58 * t71 + t59 * t70); 0.2e1 * m(7) * (t177 * t61 - t153) + 0.2e1 * m(6) * (t42 * t61 - t153) + 0.2e1 * m(5) * (-t149 + t150); m(7) * (qJD(6) * t116 + t100 * t3 + t2 * t97) + t44 * t132 + t97 * t16 - t45 * t133 + t100 * t15 + t68 - (-m(6) * pkin(4) - mrSges(5,1)) * t71 + t144; 0; 0; t103 + t128 * t88 + t104 * t87 + (m(6) * (t10 * t98 - t101 * t11) + (-t101 * t42 - t98 * t43) * mrSges(6,3) + ((m(6) * t22 + m(7) * t116 + mrSges(6,3) * t113 + t114) * t101 + (t61 * mrSges(6,3) + t35 + (m(6) + m(7)) * t21) * t98) * qJD(5)) * pkin(4) - Ifges(5,6) * t71 + Ifges(5,5) * t70 - t47 * mrSges(5,2) + t48 * mrSges(5,1); m(7) * (t177 * t87 + t43 * t88) + 0.2e1 * ((-t101 * t43 + t42 * t98) * t167 + (m(7) * (t171 * t61 - t152) / 0.2e1 + (t101 * t61 - t152) * t167) * qJD(5)) * pkin(4) + t111 - t125; 0; 0.2e1 * t88 * t79 + (-0.2e1 * t137 - 0.2e1 * t147 + 0.2e1 * t146 + (t171 * t87 + t88 * t98) * t170 + 0.2e1 * t122) * t141 + t107; -pkin(5) * t128 + pkin(9) * t104 + t103; m(7) * (pkin(9) * t177 - t165) + t111; 0; (-pkin(5) + t88) * t79 + (-t137 - t147 + m(7) * (-pkin(5) * t98 + pkin(9) * t171) + t146 + t122) * t141 + t107; -0.2e1 * pkin(5) * t79 + t107; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t130 - Ifges(7,6) * t174 + t145; t12; -t79; t89 - t120 * t101 * t141 + (t82 * t87 - t160) * qJD(6); t89 + (pkin(9) * t82 - t160) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
