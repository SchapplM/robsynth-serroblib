% Calculate joint inertia matrix for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:19
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:19:25
% EndTime: 2018-11-23 15:19:26
% DurationCPUTime: 1.23s
% Computational Cost: add. (1433->306), mult. (3431->431), div. (0->0), fcn. (3535->12), ass. (0->124)
t161 = Ifges(5,1) + Ifges(4,3);
t112 = (-pkin(3) - pkin(10));
t160 = -2 * t112;
t104 = sin(qJ(6));
t108 = cos(qJ(6));
t100 = sin(pkin(7));
t106 = sin(qJ(3));
t135 = t100 * t106;
t102 = cos(pkin(7));
t105 = sin(qJ(5));
t109 = cos(qJ(5));
t110 = cos(qJ(3));
t134 = t100 * t110;
t53 = t102 * t109 - t105 * t134;
t33 = -t104 * t53 + t108 * t135;
t52 = t102 * t105 + t109 * t134;
t17 = -mrSges(7,2) * t52 + mrSges(7,3) * t33;
t34 = t104 * t135 + t108 * t53;
t18 = mrSges(7,1) * t52 - mrSges(7,3) * t34;
t158 = -t104 * t18 + t108 * t17;
t157 = m(7) * pkin(11) + mrSges(7,3);
t67 = -mrSges(7,1) * t108 + mrSges(7,2) * t104;
t156 = m(7) * pkin(5) + mrSges(6,1) - t67;
t101 = sin(pkin(6));
t103 = cos(pkin(6));
t107 = sin(qJ(2));
t111 = cos(qJ(2));
t133 = t102 * t111;
t26 = -t103 * t134 + (t106 * t107 - t110 * t133) * t101;
t51 = -t100 * t101 * t111 + t102 * t103;
t14 = t105 * t51 - t109 * t26;
t155 = t14 ^ 2;
t28 = t103 * t135 - (-t106 * t133 - t107 * t110) * t101;
t24 = t28 ^ 2;
t154 = t33 / 0.2e1;
t153 = t34 / 0.2e1;
t69 = Ifges(7,5) * t104 + Ifges(7,6) * t108;
t152 = t69 / 0.2e1;
t151 = -t104 / 0.2e1;
t150 = t104 / 0.2e1;
t149 = t108 / 0.2e1;
t148 = pkin(2) * t100;
t147 = pkin(2) * t110;
t13 = -mrSges(7,1) * t33 + mrSges(7,2) * t34;
t38 = mrSges(6,1) * t135 - mrSges(6,3) * t53;
t146 = t13 - t38;
t126 = -pkin(3) - t147;
t79 = pkin(9) * t135;
t25 = pkin(4) * t135 + t79 + (-pkin(10) + t126) * t102;
t123 = -qJ(4) * t106 - pkin(2);
t36 = (t112 * t110 + t123) * t100;
t11 = t105 * t25 + t109 * t36;
t21 = mrSges(6,1) * t52 + mrSges(6,2) * t53;
t61 = -mrSges(5,1) * t134 - mrSges(5,3) * t102;
t145 = -t61 + t21;
t68 = mrSges(6,1) * t105 + mrSges(6,2) * t109;
t144 = t68 + mrSges(5,3);
t57 = t102 * t106 * pkin(2) + pkin(9) * t134;
t62 = mrSges(5,1) * t135 + t102 * mrSges(5,2);
t143 = t104 ^ 2 + t108 ^ 2;
t96 = t105 ^ 2;
t98 = t109 ^ 2;
t142 = t96 + t98;
t141 = Ifges(7,4) * t104;
t140 = Ifges(7,4) * t108;
t139 = qJ(4) * t28;
t136 = t109 * t14;
t132 = t104 * t109;
t131 = t105 * t112;
t130 = t108 * t109;
t129 = t109 * t112;
t5 = Ifges(7,5) * t34 + Ifges(7,6) * t33 + Ifges(7,3) * t52;
t127 = Ifges(6,5) * t53 - Ifges(6,6) * t52 + Ifges(6,3) * t135;
t41 = -t102 * qJ(4) - t57;
t125 = t142 * mrSges(6,3);
t35 = pkin(4) * t134 - t41;
t122 = Ifges(4,5) * t135 + Ifges(4,6) * t134 + t161 * t102;
t12 = pkin(5) * t52 - pkin(11) * t53 + t35;
t9 = pkin(11) * t135 + t11;
t1 = -t104 * t9 + t108 * t12;
t2 = t104 * t12 + t108 * t9;
t121 = -t1 * t104 + t108 * t2;
t16 = t105 * t26 + t109 * t51;
t3 = -t104 * t16 + t108 * t28;
t4 = t104 * t28 + t108 * t16;
t120 = -t104 * t3 + t108 * t4;
t119 = mrSges(7,1) * t104 + mrSges(7,2) * t108;
t10 = -t105 * t36 + t109 * t25;
t118 = t10 * t109 + t11 * t105;
t66 = pkin(5) * t105 - pkin(11) * t109 + qJ(4);
t39 = -t104 * t131 + t108 * t66;
t40 = t104 * t66 + t108 * t131;
t117 = -t104 * t39 + t108 * t40;
t63 = -mrSges(7,2) * t105 - mrSges(7,3) * t132;
t64 = mrSges(7,1) * t105 - mrSges(7,3) * t130;
t116 = -t104 * t64 + t108 * t63;
t115 = t105 * t16 - t136;
t44 = Ifges(7,5) * t130 - Ifges(7,6) * t132 + Ifges(7,3) * t105;
t113 = qJ(4) ^ 2;
t99 = t112 ^ 2;
t93 = Ifges(6,5) * t109;
t85 = t98 * t112;
t84 = t98 * t99;
t73 = Ifges(6,1) * t109 - Ifges(6,4) * t105;
t72 = Ifges(7,1) * t104 + t140;
t71 = Ifges(6,4) * t109 - Ifges(6,2) * t105;
t70 = Ifges(7,2) * t108 + t141;
t60 = -mrSges(4,2) * t102 + mrSges(4,3) * t134;
t59 = mrSges(4,1) * t102 - mrSges(4,3) * t135;
t58 = t119 * t109;
t56 = t102 * t147 - t79;
t55 = (-mrSges(4,1) * t110 + mrSges(4,2) * t106) * t100;
t54 = (mrSges(5,2) * t110 - mrSges(5,3) * t106) * t100;
t46 = Ifges(7,5) * t105 + (Ifges(7,1) * t108 - t141) * t109;
t45 = Ifges(7,6) * t105 + (-Ifges(7,2) * t104 + t140) * t109;
t43 = t102 * t126 + t79;
t42 = (-pkin(3) * t110 + t123) * t100;
t37 = -mrSges(6,2) * t135 - mrSges(6,3) * t52;
t20 = Ifges(6,1) * t53 - Ifges(6,4) * t52 + Ifges(6,5) * t135;
t19 = Ifges(6,4) * t53 - Ifges(6,2) * t52 + Ifges(6,6) * t135;
t8 = -pkin(5) * t135 - t10;
t7 = Ifges(7,1) * t34 + Ifges(7,4) * t33 + Ifges(7,5) * t52;
t6 = Ifges(7,4) * t34 + Ifges(7,2) * t33 + Ifges(7,6) * t52;
t15 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t155) + m(6) * (t16 ^ 2 + t155 + t24) + m(3) * (t103 ^ 2 + (t107 ^ 2 + t111 ^ 2) * t101 ^ 2) + (m(5) + m(4)) * (t26 ^ 2 + t51 ^ 2 + t24); t16 * t37 + t4 * t17 + t3 * t18 + (t54 + t55) * t51 + (-t59 + t62) * t26 + t146 * t14 + (mrSges(3,1) * t111 - mrSges(3,2) * t107) * t101 + (t60 + t145) * t28 + m(7) * (t1 * t3 + t14 * t8 + t2 * t4) + m(6) * (-t10 * t14 + t11 * t16 + t28 * t35) + m(5) * (t26 * t43 - t28 * t41 + t42 * t51) + m(4) * (-t51 * t148 - t26 * t56 + t28 * t57); 0.2e1 * t1 * t18 + 0.2e1 * t10 * t38 + 0.2e1 * t11 * t37 + 0.2e1 * t8 * t13 + 0.2e1 * t2 * t17 + t53 * t20 + 0.2e1 * t35 * t21 + t33 * t6 + t34 * t7 + 0.2e1 * t41 * t61 + 0.2e1 * t42 * t54 + 0.2e1 * t43 * t62 + 0.2e1 * t56 * t59 + 0.2e1 * t57 * t60 + Ifges(3,3) + (t5 - t19) * t52 + t122 * t102 + m(4) * (t56 ^ 2 + t57 ^ 2) + m(5) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2 + t35 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) + ((m(4) * t148 - 0.2e1 * t55) * pkin(2) + ((Ifges(5,3) + Ifges(4,2)) * t134 + (-(2 * Ifges(5,5)) + Ifges(4,6)) * t102) * t110 + ((-(2 * Ifges(5,4)) + Ifges(4,5)) * t102 + ((Ifges(4,1) + Ifges(5,2)) * t106 + 0.2e1 * (Ifges(4,4) + Ifges(5,6)) * t110) * t100 + t127) * t106) * t100; t14 * t58 + t3 * t64 + t4 * t63 + (-mrSges(4,1) + mrSges(5,2)) * t26 - t115 * mrSges(6,3) + (-mrSges(4,2) + t144) * t28 + m(7) * (-t129 * t14 + t3 * t39 + t4 * t40) + m(6) * (t112 * t115 + t139) + m(5) * (-pkin(3) * t26 + t139); m(6) * (qJ(4) * t35 + t112 * t118) - pkin(3) * t62 + t2 * t63 + t1 * t64 + t35 * t68 + t53 * t73 / 0.2e1 + t56 * mrSges(4,1) - t57 * mrSges(4,2) + t8 * t58 + t39 * t18 + t40 * t17 - t41 * mrSges(5,3) + t43 * mrSges(5,2) + t45 * t154 + t46 * t153 + t122 + (t20 / 0.2e1 + t7 * t149 + t6 * t151 - t10 * mrSges(6,3) - t146 * t112) * t109 + (-Ifges(5,5) * t110 + t106 * (-Ifges(6,6) * t105 + t93) / 0.2e1 - Ifges(5,4) * t106) * t100 + (-t71 / 0.2e1 + t44 / 0.2e1) * t52 + m(7) * (t1 * t39 - t129 * t8 + t2 * t40) + m(5) * (-pkin(3) * t43 - qJ(4) * t41) + t145 * qJ(4) + (-t19 / 0.2e1 + t5 / 0.2e1 + t112 * t37 - t11 * mrSges(6,3)) * t105; -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t39 * t64 + 0.2e1 * t40 * t63 + (t44 - t71) * t105 + m(7) * (t39 ^ 2 + t40 ^ 2 + t84) + m(6) * (t96 * t99 + t113 + t84) + m(5) * (pkin(3) ^ 2 + t113) + (-t104 * t45 + t108 * t46 + t160 * t58 + t73) * t109 + 0.2e1 * t144 * qJ(4) + t125 * t160 + t161; m(7) * (t105 * t120 - t136) + m(6) * t115 + m(5) * t26; -t146 * t109 + (t37 + t158) * t105 + m(7) * (t105 * t121 - t109 * t8) + m(6) * t118 + m(5) * t43 + t62; -m(5) * pkin(3) - t109 * t58 + mrSges(5,2) + t116 * t105 - t125 + m(7) * (t105 * t117 + t85) + m(6) * (t112 * t96 + t85); m(5) + m(6) * t142 + m(7) * (t143 * t96 + t98); -t16 * mrSges(6,2) + t120 * t157 - t14 * t156; t10 * mrSges(6,1) - t11 * mrSges(6,2) + t121 * mrSges(7,3) + t6 * t149 + t7 * t150 + t52 * t152 + t72 * t153 + t70 * t154 + t8 * t67 + t127 + (-m(7) * t8 - t13) * pkin(5) + (m(7) * t121 + t158) * pkin(11); t46 * t150 + t45 * t149 - pkin(5) * t58 + t93 + (m(7) * t117 + t116) * pkin(11) + (t112 * t156 + t72 * t149 + t70 * t151) * t109 + t117 * mrSges(7,3) + (-t112 * mrSges(6,2) - Ifges(6,6) + t152) * t105; t109 * t156 + (t143 * t157 - mrSges(6,2)) * t105; Ifges(6,3) + m(7) * (pkin(11) ^ 2 * t143 + pkin(5) ^ 2) - 0.2e1 * pkin(5) * t67 + t104 * t72 + t108 * t70 + 0.2e1 * t143 * pkin(11) * mrSges(7,3); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t5; mrSges(7,1) * t39 - mrSges(7,2) * t40 + t44; -t119 * t105; -pkin(11) * t119 + t69; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
