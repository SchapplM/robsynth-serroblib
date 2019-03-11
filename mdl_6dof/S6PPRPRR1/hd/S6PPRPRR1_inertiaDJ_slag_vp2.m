% Calculate time derivative of joint inertia matrix for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:42:07
% EndTime: 2019-03-08 18:42:11
% DurationCPUTime: 1.41s
% Computational Cost: add. (2173->272), mult. (6827->433), div. (0->0), fcn. (7220->14), ass. (0->136)
t76 = sin(qJ(6));
t79 = cos(qJ(6));
t119 = t76 ^ 2 + t79 ^ 2;
t105 = qJD(6) * t79;
t80 = cos(qJ(5));
t110 = qJD(5) * t80;
t77 = sin(qJ(5));
t84 = -t77 * t105 - t76 * t110;
t106 = qJD(6) * t77;
t100 = t76 * t106;
t101 = t79 * t110;
t156 = t100 - t101;
t155 = m(7) * pkin(10) + mrSges(7,3);
t147 = m(5) * pkin(3);
t68 = sin(pkin(13));
t154 = t147 * t68 - mrSges(5,2);
t73 = cos(pkin(12));
t74 = cos(pkin(7));
t124 = t73 * t74;
t70 = sin(pkin(7));
t78 = sin(qJ(3));
t126 = t70 * t78;
t69 = sin(pkin(12));
t71 = sin(pkin(6));
t75 = cos(pkin(6));
t81 = cos(qJ(3));
t29 = t71 * (t124 * t78 + t69 * t81) + t75 * t126;
t52 = -mrSges(7,1) * t79 + mrSges(7,2) * t76;
t153 = -m(7) * pkin(5) - mrSges(6,1) + t52;
t72 = cos(pkin(13));
t63 = -pkin(3) * t72 - pkin(4);
t152 = m(6) * t63 - mrSges(6,1) * t80 + mrSges(6,2) * t77 - t147 * t72 - mrSges(5,1);
t151 = 0.2e1 * m(7);
t62 = pkin(3) * t68 + pkin(9);
t150 = 0.2e1 * t62;
t149 = m(6) / 0.2e1;
t148 = m(7) / 0.2e1;
t28 = t70 * t75 * t81 + (t124 * t81 - t69 * t78) * t71;
t14 = t28 * t68 + t29 * t72;
t86 = -t70 * t71 * t73 + t74 * t75;
t10 = t80 * t14 + t77 * t86;
t116 = qJD(5) * t10;
t24 = t28 * qJD(3);
t25 = t29 * qJD(3);
t12 = t24 * t72 - t25 * t68;
t4 = t77 * t12 + t116;
t9 = t77 * t14 - t80 * t86;
t146 = t4 * t9;
t137 = Ifges(7,4) * t76;
t54 = Ifges(7,2) * t79 + t137;
t145 = -t54 / 0.2e1;
t144 = -t76 / 0.2e1;
t142 = pkin(5) * t77;
t141 = pkin(10) * t80;
t118 = qJD(5) * t9;
t3 = t80 * t12 - t118;
t140 = t3 * t80;
t139 = t4 * t77;
t138 = mrSges(7,3) * t77;
t136 = Ifges(7,4) * t79;
t135 = Ifges(7,5) * t76;
t134 = Ifges(7,6) * t76;
t133 = Ifges(7,6) * t79;
t132 = Ifges(7,6) * t80;
t11 = t24 * t68 + t72 * t25;
t13 = -t72 * t28 + t29 * t68;
t131 = t11 * t13;
t40 = (t68 * t81 + t72 * t78) * t70;
t32 = t40 * t77 - t80 * t74;
t115 = qJD(5) * t32;
t117 = qJD(3) * t70;
t125 = t72 * t81;
t37 = (-t68 * t78 + t125) * t117;
t17 = t80 * t37 - t115;
t130 = t17 * t80;
t33 = t40 * t80 + t74 * t77;
t114 = qJD(5) * t33;
t18 = t77 * t37 + t114;
t129 = t18 * t32;
t128 = t18 * t77;
t36 = qJD(3) * t40;
t39 = -t125 * t70 + t126 * t68;
t127 = t36 * t39;
t123 = t76 * t80;
t122 = t79 * t80;
t112 = qJD(5) * t77;
t120 = -Ifges(7,5) * t101 - Ifges(7,3) * t112;
t113 = qJD(5) * t76;
t111 = qJD(5) * t79;
t44 = -pkin(5) * t80 - pkin(10) * t77 + t63;
t30 = -t123 * t62 + t44 * t79;
t109 = qJD(6) * t30;
t31 = t122 * t62 + t44 * t76;
t108 = qJD(6) * t31;
t107 = qJD(6) * t76;
t103 = t62 * t112;
t98 = (2 * Ifges(6,4)) + t134;
t51 = (-t141 + t142) * qJD(5);
t15 = -t103 * t79 + t51 * t76 + t109;
t97 = t15 - t109;
t5 = -t10 * t76 + t13 * t79;
t1 = qJD(6) * t5 + t11 * t76 + t3 * t79;
t6 = t10 * t79 + t13 * t76;
t2 = -qJD(6) * t6 + t11 * t79 - t3 * t76;
t96 = t1 * t79 - t2 * t76;
t95 = t18 * t9 + t32 * t4;
t19 = -t33 * t76 + t39 * t79;
t7 = qJD(6) * t19 + t17 * t79 + t36 * t76;
t20 = t33 * t79 + t39 * t76;
t8 = -qJD(6) * t20 - t17 * t76 + t36 * t79;
t94 = t7 * t79 - t76 * t8;
t93 = t77 * mrSges(6,1) + t80 * mrSges(6,2);
t92 = mrSges(7,1) * t76 + mrSges(7,2) * t79;
t91 = Ifges(7,1) * t79 - t137;
t55 = Ifges(7,1) * t76 + t136;
t90 = -Ifges(7,2) * t76 + t136;
t89 = t11 * t39 + t13 * t36;
t88 = t110 * t9 + t139;
t85 = t110 * t32 + t128;
t65 = Ifges(7,5) * t105;
t50 = -mrSges(7,1) * t80 - t138 * t79;
t49 = mrSges(7,2) * t80 - t138 * t76;
t48 = t91 * qJD(6);
t47 = t90 * qJD(6);
t46 = t93 * qJD(5);
t45 = t92 * qJD(6);
t43 = t92 * t77;
t42 = -Ifges(7,5) * t80 + t77 * t91;
t41 = t77 * t90 - t132;
t35 = -mrSges(7,2) * t112 + mrSges(7,3) * t84;
t34 = mrSges(7,1) * t112 + t156 * mrSges(7,3);
t27 = t84 * mrSges(7,1) + t156 * mrSges(7,2);
t22 = -t55 * t106 + (Ifges(7,5) * t77 + t80 * t91) * qJD(5);
t21 = -t54 * t106 + (Ifges(7,6) * t77 + t80 * t90) * qJD(5);
t16 = t103 * t76 + t51 * t79 - t108;
t23 = [0.2e1 * m(7) * (t1 * t6 + t2 * t5 + t146) + 0.2e1 * m(6) * (t10 * t3 + t131 + t146) + 0.2e1 * m(4) * (t24 * t29 - t25 * t28) + 0.2e1 * m(5) * (t12 * t14 + t131); m(7) * (t1 * t20 + t19 * t2 + t5 * t8 + t6 * t7 + t95) + m(6) * (t10 * t17 + t3 * t33 + t89 + t95) + m(5) * (t12 * t40 + t14 * t37 + t89) + m(4) * (t24 * t78 - t25 * t81 + (-t28 * t78 + t29 * t81) * qJD(3)) * t70; 0.2e1 * m(7) * (t19 * t8 + t20 * t7 + t129) + 0.2e1 * m(6) * (t17 * t33 + t127 + t129) + 0.2e1 * m(5) * (t37 * t40 + t127); -t25 * mrSges(4,1) - t24 * mrSges(4,2) + t1 * t49 + t13 * t46 + t2 * t50 - t9 * t27 + t5 * t34 + t6 * t35 + t4 * t43 + m(7) * (t1 * t31 + t15 * t6 + t16 * t5 + t2 * t30) + (t88 * t148 + (-t10 * t112 + t140 + t88) * t149) * t150 + (t140 + t139 + (-t10 * t77 + t80 * t9) * qJD(5)) * mrSges(6,3) + t154 * t12 + t152 * t11; t18 * t43 + t19 * t34 + t20 * t35 - t32 * t27 + t39 * t46 + t7 * t49 + t8 * t50 + (-mrSges(4,1) * t78 - mrSges(4,2) * t81) * t117 + m(7) * (t15 * t20 + t16 * t19 + t30 * t8 + t31 * t7) + (t85 * t148 + (-t112 * t33 + t130 + t85) * t149) * t150 + (t130 + t128 + (t32 * t80 - t33 * t77) * qJD(5)) * mrSges(6,3) + t154 * t37 + t152 * t36; (t15 * t31 + t16 * t30) * t151 + 0.2e1 * t15 * t49 + 0.2e1 * t31 * t35 + 0.2e1 * t16 * t50 + 0.2e1 * t30 * t34 + 0.2e1 * t63 * t46 + ((t150 * t43 - t76 * t41 + t79 * t42 + t80 * t98) * qJD(5) + t120) * t80 + (-t76 * t21 + t79 * t22 - 0.2e1 * t62 * t27 + (-t80 * (-t133 - t135) - t79 * t41 - t76 * t42) * qJD(6) + ((Ifges(7,5) * t79 - t98) * t77 + (t151 * t62 ^ 2 + (2 * Ifges(6,1)) - (2 * Ifges(6,2)) - Ifges(7,3)) * t80) * qJD(5)) * t77; 0.2e1 * ((t111 * t6 - t113 * t5 - t4) * t148 + (-t4 + t116) * t149) * t80 + 0.2e1 * ((-t105 * t5 - t107 * t6 + t118 + t96) * t148 + (t3 + t118) * t149) * t77; 0.2e1 * ((t111 * t20 - t113 * t19 - t18) * t148 + (-t18 + t114) * t149) * t80 + 0.2e1 * ((-t105 * t19 - t107 * t20 + t115 + t94) * t148 + (t17 + t115) * t149) * t77; t80 * t27 + (m(7) * (-t105 * t30 - t107 * t31 + t15 * t79 - t16 * t76) - t49 * t107 + t79 * t35 - t50 * t105 - t76 * t34) * t77 + (m(7) * (t31 * t122 - t30 * t123 + (t77 ^ 2 - t80 ^ 2) * t62) + t49 * t122 - t50 * t123 + t77 * t43) * qJD(5); (-0.1e1 + t119) * t77 * t110 * t151; -t3 * mrSges(6,2) + t9 * t45 + t155 * ((-t5 * t79 - t6 * t76) * qJD(6) + t96) + t153 * t4; -t17 * mrSges(6,2) + t32 * t45 + t155 * ((-t19 * t79 - t20 * t76) * qJD(6) + t94) + t153 * t18; pkin(5) * t27 + (-t65 / 0.2e1 + (t153 * t62 + Ifges(6,5)) * qJD(5)) * t80 + (t110 * t145 - t16 * mrSges(7,3) + t22 / 0.2e1 + (-t31 * mrSges(7,3) - t41 / 0.2e1 + t132 / 0.2e1) * qJD(6) + (m(7) * (-t16 - t108) - qJD(6) * t49 - t34) * pkin(10)) * t76 + (t55 * t110 / 0.2e1 + qJD(6) * t42 / 0.2e1 + t21 / 0.2e1 + t97 * mrSges(7,3) + (m(7) * t97 - qJD(6) * t50 + t35) * pkin(10)) * t79 + (t62 * t45 + t47 * t144 + t79 * t48 / 0.2e1 + (t144 * t55 + t145 * t79) * qJD(6) + (t62 * mrSges(6,2) - Ifges(6,6) + t135 / 0.2e1 + t133 / 0.2e1) * qJD(5)) * t77; -t80 * t45 + (t77 * t52 + m(7) * (t119 * t141 - t142) + t119 * t80 * mrSges(7,3) - t93) * qJD(5); -0.2e1 * pkin(5) * t45 + t47 * t79 + t48 * t76 + (-t54 * t76 + t55 * t79) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1; mrSges(7,1) * t8 - mrSges(7,2) * t7; mrSges(7,1) * t16 - mrSges(7,2) * t15 - Ifges(7,5) * t100 + Ifges(7,6) * t84 - t120; t27; t65 + (pkin(10) * t52 - t134) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t23(1) t23(2) t23(4) t23(7) t23(11) t23(16); t23(2) t23(3) t23(5) t23(8) t23(12) t23(17); t23(4) t23(5) t23(6) t23(9) t23(13) t23(18); t23(7) t23(8) t23(9) t23(10) t23(14) t23(19); t23(11) t23(12) t23(13) t23(14) t23(15) t23(20); t23(16) t23(17) t23(18) t23(19) t23(20) t23(21);];
Mq  = res;
