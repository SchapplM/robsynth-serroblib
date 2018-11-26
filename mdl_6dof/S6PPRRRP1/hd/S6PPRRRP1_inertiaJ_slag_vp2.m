% Calculate joint inertia matrix for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2018-11-23 14:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:50:45
% EndTime: 2018-11-23 14:50:46
% DurationCPUTime: 0.89s
% Computational Cost: add. (816->220), mult. (2130->311), div. (0->0), fcn. (2267->12), ass. (0->97)
t131 = Ifges(6,5) + Ifges(7,5);
t99 = Ifges(6,6) + Ifges(7,6);
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t39 = -t77 * mrSges(7,1) + t74 * mrSges(7,2);
t55 = -t77 * pkin(5) - pkin(4);
t130 = m(7) * t55 + t39;
t40 = -t77 * mrSges(6,1) + t74 * mrSges(6,2);
t129 = -m(6) * pkin(4) - mrSges(5,1) + t40;
t128 = 2 * pkin(9);
t127 = -2 * mrSges(7,3);
t126 = m(6) + m(7);
t123 = t131 * t74 + t99 * t77;
t122 = m(6) * pkin(10) + mrSges(6,3) + mrSges(7,3);
t121 = t129 + t130;
t71 = cos(pkin(12));
t72 = cos(pkin(7));
t105 = t71 * t72;
t69 = sin(pkin(7));
t76 = sin(qJ(3));
t107 = t69 * t76;
t68 = sin(pkin(12));
t70 = sin(pkin(6));
t73 = cos(pkin(6));
t79 = cos(qJ(3));
t12 = t73 * t107 + (t76 * t105 + t68 * t79) * t70;
t23 = -t69 * t70 * t71 + t73 * t72;
t75 = sin(qJ(4));
t78 = cos(qJ(4));
t6 = t75 * t12 - t78 * t23;
t120 = t6 ^ 2;
t106 = t69 * t79;
t10 = -t73 * t106 + (-t105 * t79 + t68 * t76) * t70;
t119 = t10 ^ 2;
t24 = t75 * t107 - t78 * t72;
t118 = t24 ^ 2;
t117 = m(7) * pkin(5);
t115 = pkin(9) * t78;
t114 = t24 * t6;
t113 = t6 * t75;
t112 = Ifges(6,4) * t74;
t111 = Ifges(6,4) * t77;
t110 = Ifges(7,4) * t74;
t109 = Ifges(7,4) * t77;
t108 = t24 * t75;
t104 = t74 * t75;
t103 = t75 * t77;
t102 = t78 * mrSges(5,3);
t41 = -t78 * mrSges(5,1) + t75 * mrSges(5,2);
t101 = mrSges(4,1) - t41;
t100 = -mrSges(6,2) - mrSges(7,2);
t98 = Ifges(6,3) + Ifges(7,3);
t97 = -qJ(6) - pkin(10);
t32 = t78 * mrSges(7,2) - mrSges(7,3) * t104;
t33 = t78 * mrSges(6,2) - mrSges(6,3) * t104;
t96 = t32 + t33;
t34 = -t78 * mrSges(7,1) - mrSges(7,3) * t103;
t35 = -t78 * mrSges(6,1) - mrSges(6,3) * t103;
t95 = t34 + t35;
t27 = mrSges(7,1) * t104 + mrSges(7,2) * t103;
t93 = t131 * t103;
t37 = -t78 * pkin(4) - t75 * pkin(10) - pkin(3);
t17 = t77 * t115 + t74 * t37;
t90 = t74 ^ 2 + t77 ^ 2;
t89 = qJ(6) * t75;
t86 = mrSges(6,1) + mrSges(7,1) + t117;
t83 = mrSges(6,1) * t74 + mrSges(6,2) * t77;
t28 = t83 * t75;
t85 = t75 * mrSges(5,3) + t27 + t28;
t81 = pkin(9) ^ 2;
t67 = t78 ^ 2;
t65 = t75 ^ 2;
t62 = t69 ^ 2;
t61 = t65 * t81;
t53 = t62 * t79 ^ 2;
t46 = Ifges(6,1) * t74 + t111;
t45 = Ifges(7,1) * t74 + t109;
t44 = Ifges(6,2) * t77 + t112;
t43 = Ifges(7,2) * t77 + t110;
t42 = t97 * t77;
t38 = t97 * t74;
t36 = (pkin(5) * t74 + pkin(9)) * t75;
t30 = t77 * t37;
t26 = t78 * t107 + t75 * t72;
t22 = -Ifges(6,5) * t78 + (Ifges(6,1) * t77 - t112) * t75;
t21 = -Ifges(7,5) * t78 + (Ifges(7,1) * t77 - t110) * t75;
t20 = -Ifges(6,6) * t78 + (-Ifges(6,2) * t74 + t111) * t75;
t19 = -Ifges(7,6) * t78 + (-Ifges(7,2) * t74 + t109) * t75;
t16 = -t74 * t115 + t30;
t15 = -t74 * t89 + t17;
t14 = -t74 * t106 + t77 * t26;
t13 = -t77 * t106 - t74 * t26;
t9 = -t77 * t89 + t30 + (-pkin(9) * t74 - pkin(5)) * t78;
t8 = t78 * t12 + t75 * t23;
t4 = t74 * t10 + t77 * t8;
t3 = t77 * t10 - t74 * t8;
t1 = [m(2) + m(5) * (t8 ^ 2 + t119 + t120) + m(4) * (t12 ^ 2 + t23 ^ 2 + t119) + m(3) * (t73 ^ 2 + (t68 ^ 2 + t71 ^ 2) * t70 ^ 2) + (t3 ^ 2 + t4 ^ 2 + t120) * t126; m(3) * t73 + m(5) * (-t10 * t106 + t26 * t8 + t114) + m(4) * (t72 * t23 + (-t10 * t79 + t12 * t76) * t69) + (t13 * t3 + t14 * t4 + t114) * t126; m(3) + m(5) * (t26 ^ 2 + t118 + t53) + m(4) * (t62 * t76 ^ 2 + t72 ^ 2 + t53) + (t13 ^ 2 + t14 ^ 2 + t118) * t126; t8 * t102 - t12 * mrSges(4,2) + t96 * t4 + t95 * t3 - t101 * t10 + t85 * t6 + m(6) * (pkin(9) * t113 + t16 * t3 + t17 * t4) + m(7) * (t15 * t4 + t9 * t3 + t36 * t6) + m(5) * (-pkin(3) * t10 + (t8 * t78 + t113) * pkin(9)); t26 * t102 + t96 * t14 + t95 * t13 + (-t76 * mrSges(4,2) + t101 * t79) * t69 + t85 * t24 + m(7) * (t9 * t13 + t15 * t14 + t36 * t24) + m(6) * (pkin(9) * t108 + t16 * t13 + t17 * t14) + m(5) * (pkin(3) * t106 + (t26 * t78 + t108) * pkin(9)); -0.2e1 * pkin(3) * t41 + 0.2e1 * t15 * t32 + 0.2e1 * t16 * t35 + 0.2e1 * t17 * t33 + 0.2e1 * t36 * t27 + 0.2e1 * t9 * t34 + Ifges(4,3) + (t65 + t67) * mrSges(5,3) * t128 + m(7) * (t15 ^ 2 + t36 ^ 2 + t9 ^ 2) + m(6) * (t16 ^ 2 + t17 ^ 2 + t61) + m(5) * (pkin(3) ^ 2 + t67 * t81 + t61) + ((Ifges(5,2) + t98) * t78 - t93) * t78 + (Ifges(5,1) * t75 + 0.2e1 * Ifges(5,4) * t78 + t28 * t128 + (t21 + t22) * t77 + (t99 * t78 - t19 - t20) * t74) * t75; -t8 * mrSges(5,2) + m(7) * (t38 * t3 - t42 * t4) + t121 * t6 + t122 * (-t3 * t74 + t4 * t77); -t26 * mrSges(5,2) + m(7) * (t38 * t13 - t42 * t14) + t121 * t24 + t122 * (-t13 * t74 + t14 * t77); t55 * t27 + t38 * t34 + t36 * t39 - t42 * t32 - pkin(4) * t28 + m(7) * (-t42 * t15 + t55 * t36 + t38 * t9) + (t17 * mrSges(6,3) + t15 * mrSges(7,3) + t19 / 0.2e1 + t20 / 0.2e1 + (m(6) * t17 + t33) * pkin(10)) * t77 + (-t16 * mrSges(6,3) - t9 * mrSges(7,3) + t22 / 0.2e1 + t21 / 0.2e1 + (-m(6) * t16 - t35) * pkin(10)) * t74 + (Ifges(5,5) + t129 * pkin(9) + (t45 / 0.2e1 + t46 / 0.2e1) * t77 + (-t43 / 0.2e1 - t44 / 0.2e1) * t74) * t75 + (Ifges(5,6) - pkin(9) * mrSges(5,2) - t123 / 0.2e1) * t78; -0.2e1 * pkin(4) * t40 + 0.2e1 * t55 * t39 + Ifges(5,3) + 0.2e1 * t90 * pkin(10) * mrSges(6,3) + m(7) * (t38 ^ 2 + t42 ^ 2 + t55 ^ 2) + m(6) * (t90 * pkin(10) ^ 2 + pkin(4) ^ 2) + (t42 * t127 + t43 + t44) * t77 + (t38 * t127 + t45 + t46) * t74; t100 * t4 + t86 * t3; t100 * t14 + t86 * t13; t16 * mrSges(6,1) + t9 * mrSges(7,1) - t17 * mrSges(6,2) - t15 * mrSges(7,2) - t98 * t78 - t99 * t104 + (m(7) * t9 + t34) * pkin(5) + t93; t38 * mrSges(7,1) + t42 * mrSges(7,2) - t83 * pkin(10) + (m(7) * t38 - t74 * mrSges(7,3)) * pkin(5) + t123; (0.2e1 * mrSges(7,1) + t117) * pkin(5) + t98; m(7) * t6; m(7) * t24; m(7) * t36 + t27; t130; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
