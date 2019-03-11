% Calculate time derivative of joint inertia matrix for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:06
% EndTime: 2019-03-09 02:50:09
% DurationCPUTime: 1.57s
% Computational Cost: add. (2983->273), mult. (6582->404), div. (0->0), fcn. (6360->8), ass. (0->120)
t156 = 2 * mrSges(5,1) + 2 * mrSges(4,3);
t101 = sin(pkin(10));
t103 = cos(pkin(10));
t106 = sin(qJ(6));
t108 = cos(qJ(6));
t114 = t101 * t106 - t103 * t108;
t115 = t101 * t108 + t103 * t106;
t75 = t114 * qJD(6);
t76 = t115 * qJD(6);
t155 = t114 * t76 - t115 * t75;
t99 = t103 ^ 2;
t124 = (t101 ^ 2 + t99) * qJD(5);
t102 = sin(pkin(9));
t107 = sin(qJ(3));
t151 = cos(qJ(3));
t125 = qJD(3) * t151;
t104 = cos(pkin(9));
t126 = t151 * t104;
t144 = pkin(7) + qJ(2);
t87 = t144 * t102;
t89 = t144 * t104;
t45 = (qJD(2) * t102 + qJD(3) * t89) * t107 - qJD(2) * t126 + t87 * t125;
t154 = 2 * m(6);
t153 = 2 * m(7);
t133 = t102 * t107;
t77 = qJD(3) * t133 - t104 * t125;
t84 = t102 * t151 + t104 * t107;
t116 = qJ(4) * t77 - qJD(4) * t84;
t78 = t84 * qJD(3);
t39 = pkin(3) * t78 + t116;
t152 = -0.2e1 * t39;
t150 = Ifges(6,5) * t77;
t149 = Ifges(6,6) * t77;
t146 = mrSges(5,2) - mrSges(4,1);
t145 = pkin(3) + qJ(5);
t143 = -pkin(8) - t145;
t82 = -t126 + t133;
t26 = qJD(5) * t82 + t145 * t78 + t116;
t68 = -t107 * t87 + t151 * t89;
t46 = qJD(2) * t84 + qJD(3) * t68;
t32 = -t77 * pkin(4) + t46;
t11 = t101 * t32 + t103 * t26;
t128 = -pkin(2) * t104 - pkin(1);
t112 = -qJ(4) * t84 + t128;
t40 = t145 * t82 + t112;
t67 = t107 * t89 + t151 * t87;
t47 = pkin(4) * t84 + t67;
t17 = t101 * t47 + t103 * t40;
t142 = -Ifges(7,5) * t76 + Ifges(7,6) * t75;
t139 = Ifges(6,4) * t101;
t138 = Ifges(6,4) * t103;
t137 = t101 * Ifges(6,1);
t136 = t101 * t78;
t135 = t103 * t78;
t134 = t103 * t82;
t130 = 2 * mrSges(7,3);
t50 = t114 * t82;
t24 = -qJD(6) * t50 + t115 * t78;
t25 = -t114 * t78 - t76 * t82;
t129 = Ifges(7,5) * t24 + Ifges(7,6) * t25 - Ifges(7,3) * t77;
t127 = -pkin(5) * t103 - pkin(4);
t8 = -t25 * mrSges(7,1) + mrSges(7,2) * t24;
t49 = -mrSges(6,1) * t135 + mrSges(6,2) * t136;
t55 = -mrSges(7,1) * t75 - mrSges(7,2) * t76;
t120 = -t45 * t68 + t46 * t67;
t119 = -Ifges(6,5) * t101 - Ifges(6,6) * t103;
t30 = t103 * t32;
t10 = -t101 * t26 + t30;
t118 = t10 * t103 + t101 * t11;
t44 = t103 * t47;
t12 = pkin(5) * t84 + t44 + (-pkin(8) * t82 - t40) * t101;
t13 = pkin(8) * t134 + t17;
t3 = -t106 * t13 + t108 * t12;
t4 = t106 * t12 + t108 * t13;
t85 = t143 * t101;
t86 = t143 * t103;
t62 = t106 * t86 + t108 * t85;
t61 = -t106 * t85 + t108 * t86;
t117 = t155 * t153;
t7 = -pkin(5) * t77 + t30 + (-pkin(8) * t78 - t26) * t101;
t9 = pkin(8) * t135 + t11;
t1 = qJD(6) * t3 + t106 * t7 + t108 * t9;
t2 = -qJD(6) * t4 - t106 * t9 + t108 * t7;
t111 = -t1 * t115 + t114 * t2 + t3 * t76 + t4 * t75;
t41 = -qJD(5) * t115 + qJD(6) * t61;
t42 = qJD(5) * t114 - qJD(6) * t62;
t109 = t114 * t42 - t115 * t41 + t61 * t76 + t62 * t75;
t93 = pkin(5) * t101 + qJ(4);
t88 = mrSges(6,1) * t101 + mrSges(6,2) * t103;
t71 = t77 * mrSges(4,2);
t70 = t77 * mrSges(5,3);
t65 = -Ifges(7,1) * t114 - Ifges(7,4) * t115;
t64 = -Ifges(7,4) * t114 - Ifges(7,2) * t115;
t63 = mrSges(7,1) * t115 - mrSges(7,2) * t114;
t60 = -mrSges(6,2) * t84 + mrSges(6,3) * t134;
t59 = -mrSges(6,3) * t101 * t82 + mrSges(6,1) * t84;
t58 = pkin(3) * t82 + t112;
t57 = -Ifges(7,1) * t76 + Ifges(7,4) * t75;
t56 = -Ifges(7,4) * t76 + Ifges(7,2) * t75;
t54 = (-mrSges(6,1) * t103 + mrSges(6,2) * t101) * t82;
t53 = mrSges(6,2) * t77 + mrSges(6,3) * t135;
t52 = -mrSges(6,1) * t77 - mrSges(6,3) * t136;
t51 = t115 * t82;
t48 = -pkin(4) * t82 + t68;
t37 = mrSges(7,1) * t84 - mrSges(7,3) * t51;
t36 = -mrSges(7,2) * t84 - mrSges(7,3) * t50;
t35 = -t150 + (t137 + t138) * t78;
t34 = -t149 + (t103 * Ifges(6,2) + t139) * t78;
t33 = t127 * t82 + t68;
t31 = -pkin(4) * t78 - t45;
t28 = mrSges(7,1) * t50 + mrSges(7,2) * t51;
t27 = t127 * t78 - t45;
t19 = Ifges(7,1) * t51 - Ifges(7,4) * t50 + Ifges(7,5) * t84;
t18 = Ifges(7,4) * t51 - Ifges(7,2) * t50 + Ifges(7,6) * t84;
t16 = -t101 * t40 + t44;
t15 = mrSges(7,2) * t77 + mrSges(7,3) * t25;
t14 = -mrSges(7,1) * t77 - mrSges(7,3) * t24;
t6 = Ifges(7,1) * t24 + Ifges(7,4) * t25 - Ifges(7,5) * t77;
t5 = Ifges(7,4) * t24 + Ifges(7,2) * t25 - Ifges(7,6) * t77;
t20 = [(mrSges(5,3) * t152 + t156 * t46 + t129) * t84 + (mrSges(5,2) * t152 + t101 * t35 + t103 * t34 + t156 * t45) * t82 + (-0.2e1 * t58 * mrSges(5,2) + 0.2e1 * t128 * mrSges(4,1) - t68 * t156 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6) - t119) * t84 + (t99 * Ifges(6,2) + (2 * Ifges(4,2)) + (2 * Ifges(5,3)) + (t137 + 0.2e1 * t138) * t101) * t82) * t78 + (-Ifges(7,5) * t51 + Ifges(7,6) * t50 - t67 * t156 + (-(2 * Ifges(4,1)) - (2 * Ifges(5,2)) - (2 * Ifges(6,3)) - Ifges(7,3)) * t84 + (0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + t119) * t82) * t77 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t102 ^ 2 + t104 ^ 2) * qJD(2) + (t1 * t4 + t2 * t3 + t27 * t33) * t153 + (t10 * t16 + t11 * t17 + t31 * t48) * t154 - 0.2e1 * t128 * t71 + 0.2e1 * m(5) * (t39 * t58 + t120) + 0.2e1 * m(4) * t120 + 0.2e1 * t3 * t14 + 0.2e1 * t4 * t15 + t24 * t19 + t25 * t18 + 0.2e1 * t27 * t28 + 0.2e1 * t33 * t8 + 0.2e1 * t1 * t36 + 0.2e1 * t2 * t37 + 0.2e1 * t48 * t49 - t50 * t5 + t51 * t6 + 0.2e1 * t16 * t52 + 0.2e1 * t17 * t53 + 0.2e1 * t31 * t54 + 0.2e1 * t10 * t59 + 0.2e1 * t11 * t60 + 0.2e1 * t58 * t70; -t101 * t52 + t103 * t53 - t115 * t14 - t114 * t15 - t76 * t36 + t75 * t37 + t70 - t71 - t146 * t78 + m(7) * (-t1 * t114 - t115 * t2 + t3 * t75 - t4 * t76) + m(6) * (-t10 * t101 + t103 * t11) + m(5) * t39; t117; -t115 * t5 / 0.2e1 + (pkin(3) * mrSges(5,1) + Ifges(5,4) - Ifges(4,5) + Ifges(7,5) * t114 / 0.2e1 + Ifges(7,6) * t115 / 0.2e1) * t77 - t114 * t6 / 0.2e1 + (-t11 * mrSges(6,3) - qJD(5) * t60 - t145 * t53 - t34 / 0.2e1 + t149 / 0.2e1) * t101 + (-t10 * mrSges(6,3) - qJD(5) * t59 - t145 * t52 + t35 / 0.2e1 - t150 / 0.2e1) * t103 + m(6) * (qJ(4) * t31 + qJD(4) * t48 - t118 * t145 + (-t101 * t17 - t103 * t16) * qJD(5)) + m(5) * (-pkin(3) * t46 - qJ(4) * t45 + qJD(4) * t68) + m(7) * (qJD(4) * t33 + t1 * t62 + t2 * t61 + t27 * t93 + t3 * t42 + t4 * t41) + t146 * t46 + t84 * t142 / 0.2e1 + (-qJ(4) * mrSges(5,1) + t101 * (Ifges(6,1) * t103 - t139) / 0.2e1 + t103 * (-Ifges(6,2) * t101 + t138) / 0.2e1 + Ifges(5,5) - Ifges(4,6)) * t78 + t111 * mrSges(7,3) + t41 * t36 + t42 * t37 + qJ(4) * t49 + t33 * t55 - t50 * t56 / 0.2e1 + t51 * t57 / 0.2e1 + t61 * t14 + t62 * t15 + t27 * t63 + t25 * t64 / 0.2e1 + t24 * t65 / 0.2e1 + t75 * t18 / 0.2e1 - t76 * t19 / 0.2e1 + t31 * t88 + t93 * t8 + (-mrSges(5,3) + mrSges(4,2)) * t45 + (-t82 * mrSges(5,1) + t28 + t54) * qJD(4); m(7) * (-t114 * t41 - t115 * t42 + t61 * t75 - t62 * t76); 0.2e1 * t93 * t55 - t115 * t56 - t114 * t57 + t75 * t64 - t76 * t65 + (qJD(4) * t93 + t41 * t62 + t42 * t61) * t153 + (qJ(4) * qJD(4) + t124 * t145) * t154 + t109 * t130 + 0.2e1 * mrSges(6,3) * t124 + 0.2e1 * (m(5) * qJ(4) + mrSges(5,3) + t63 + t88) * qJD(4); m(5) * t46 + m(6) * t118 - m(7) * t111 - t77 * mrSges(5,1) + t101 * t53 + t103 * t52 - t114 * t14 + t115 * t15 - t75 * t36 - t76 * t37; 0; -m(6) * t124 - m(7) * t109 - t130 * t155; t117; m(6) * t31 + m(7) * t27 + t49 + t8; 0; (m(6) + m(7)) * qJD(4) + t55; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t129; -t55; mrSges(7,1) * t42 - mrSges(7,2) * t41 + t142; -mrSges(7,1) * t76 + mrSges(7,2) * t75; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
