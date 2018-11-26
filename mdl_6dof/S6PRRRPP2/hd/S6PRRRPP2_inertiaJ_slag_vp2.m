% Calculate joint inertia matrix for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2018-11-23 15:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:21:02
% EndTime: 2018-11-23 15:21:03
% DurationCPUTime: 0.97s
% Computational Cost: add. (695->272), mult. (1530->341), div. (0->0), fcn. (1290->8), ass. (0->103)
t133 = m(6) + m(7);
t79 = sin(qJ(4));
t82 = cos(qJ(4));
t132 = t79 ^ 2 + t82 ^ 2;
t77 = sin(pkin(6));
t84 = cos(qJ(2));
t111 = t77 * t84;
t81 = sin(qJ(2));
t112 = t77 * t81;
t78 = cos(pkin(6));
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t26 = t112 * t83 + t78 * t80;
t7 = t111 * t82 + t26 * t79;
t9 = -t111 * t79 + t26 * t82;
t131 = t7 * t79 + t82 * t9;
t130 = 2 * pkin(8);
t24 = t112 * t80 - t78 * t83;
t128 = t24 ^ 2;
t127 = -2 * mrSges(7,3);
t126 = pkin(4) + pkin(5);
t125 = pkin(8) * t83;
t121 = Ifges(5,4) * t79;
t120 = Ifges(5,4) * t82;
t119 = Ifges(7,4) * t79;
t118 = Ifges(7,4) * t82;
t117 = Ifges(6,5) * t79;
t116 = Ifges(6,5) * t82;
t115 = Ifges(7,5) * t83;
t114 = t24 * t80;
t113 = t26 * t83;
t110 = t79 * t80;
t109 = t80 * t82;
t108 = -mrSges(6,1) - mrSges(7,1);
t107 = mrSges(6,2) - mrSges(7,3);
t106 = mrSges(6,3) + mrSges(7,2);
t105 = Ifges(5,6) + Ifges(7,6);
t104 = pkin(9) - qJ(6);
t33 = mrSges(5,2) * t83 - mrSges(5,3) * t110;
t37 = -mrSges(6,2) * t110 - mrSges(6,3) * t83;
t103 = t33 + t37;
t35 = -mrSges(5,1) * t83 - mrSges(5,3) * t109;
t36 = t83 * mrSges(6,1) + mrSges(6,2) * t109;
t102 = -t35 + t36;
t43 = -mrSges(5,1) * t82 + mrSges(5,2) * t79;
t101 = t43 - mrSges(4,1);
t39 = -pkin(3) * t83 - pkin(9) * t80 - pkin(2);
t14 = t82 * t125 + t79 * t39;
t42 = t82 * mrSges(7,1) + t79 * mrSges(7,2);
t100 = t132 * pkin(9) ^ 2;
t99 = qJ(5) * t82;
t98 = qJ(6) * t80;
t96 = Ifges(5,3) + Ifges(7,3) + Ifges(6,2);
t95 = -Ifges(6,6) * t110 + (-Ifges(6,4) - Ifges(5,5)) * t109;
t94 = qJ(5) * t79 + pkin(3);
t57 = t79 * t125;
t13 = t39 * t82 - t57;
t93 = t131 * pkin(9);
t28 = -mrSges(7,1) * t110 + mrSges(7,2) * t109;
t34 = t83 * mrSges(7,1) - mrSges(7,3) * t109;
t11 = -qJ(5) * t83 + t14;
t91 = t79 * mrSges(5,1) + t82 * mrSges(5,2);
t90 = t79 * mrSges(6,1) - t82 * mrSges(6,3);
t89 = -pkin(4) * t79 + t99;
t88 = pkin(8) ^ 2;
t86 = qJ(5) ^ 2;
t76 = t83 ^ 2;
t74 = t80 ^ 2;
t72 = t77 ^ 2;
t71 = t83 * pkin(4);
t69 = t74 * t88;
t67 = Ifges(6,4) * t79;
t66 = Ifges(5,5) * t79;
t65 = Ifges(5,6) * t82;
t59 = t72 * t84 ^ 2;
t51 = Ifges(5,1) * t79 + t120;
t50 = Ifges(6,1) * t79 - t116;
t49 = Ifges(7,1) * t79 - t118;
t48 = Ifges(5,2) * t82 + t121;
t47 = -Ifges(7,2) * t82 + t119;
t46 = -Ifges(6,3) * t82 + t117;
t45 = t104 * t82;
t44 = -mrSges(4,1) * t83 + mrSges(4,2) * t80;
t41 = -mrSges(6,1) * t82 - mrSges(6,3) * t79;
t40 = t104 * t79;
t38 = -pkin(4) * t82 - t94;
t32 = -mrSges(7,2) * t83 + mrSges(7,3) * t110;
t30 = t126 * t82 + t94;
t29 = t91 * t80;
t27 = t90 * t80;
t21 = (pkin(8) - t89) * t80;
t20 = -Ifges(5,5) * t83 + (Ifges(5,1) * t82 - t121) * t80;
t19 = -Ifges(6,4) * t83 + (Ifges(6,1) * t82 + t117) * t80;
t18 = t115 + (Ifges(7,1) * t82 + t119) * t80;
t17 = -Ifges(5,6) * t83 + (-Ifges(5,2) * t79 + t120) * t80;
t16 = Ifges(7,6) * t83 + (Ifges(7,2) * t79 + t118) * t80;
t15 = -Ifges(6,6) * t83 + (Ifges(6,3) * t79 + t116) * t80;
t12 = -t13 + t71;
t10 = (-t126 * t79 - pkin(8) + t99) * t80;
t4 = t9 * qJ(5);
t3 = t79 * t98 + t11;
t1 = pkin(5) * t83 + t57 + t71 + (-t39 - t98) * t82;
t2 = [m(2) + m(4) * (t26 ^ 2 + t128 + t59) + m(3) * (t72 * t81 ^ 2 + t78 ^ 2 + t59) + (m(5) + t133) * (t7 ^ 2 + t9 ^ 2 + t128); mrSges(4,3) * t113 + (t32 + t103) * t9 + (-t81 * mrSges(3,2) + (mrSges(3,1) - t44) * t84) * t77 + (t34 + t102) * t7 + (t80 * mrSges(4,3) + t27 - t28 + t29) * t24 + m(6) * (t11 * t9 + t12 * t7 + t21 * t24) + m(7) * (t1 * t7 - t10 * t24 + t3 * t9) + m(5) * (pkin(8) * t114 - t13 * t7 + t14 * t9) + m(4) * (pkin(2) * t111 + (t113 + t114) * pkin(8)); -0.2e1 * pkin(2) * t44 + 0.2e1 * t1 * t34 + 0.2e1 * t10 * t28 + 0.2e1 * t11 * t37 + 0.2e1 * t12 * t36 + 0.2e1 * t13 * t35 + 0.2e1 * t14 * t33 + 0.2e1 * t21 * t27 + 0.2e1 * t3 * t32 + Ifges(3,3) + (t74 + t76) * mrSges(4,3) * t130 + m(4) * (pkin(2) ^ 2 + t76 * t88 + t69) + m(5) * (t13 ^ 2 + t14 ^ 2 + t69) + m(6) * (t11 ^ 2 + t12 ^ 2 + t21 ^ 2) + m(7) * (t1 ^ 2 + t10 ^ 2 + t3 ^ 2) + ((Ifges(4,2) + t96) * t83 + t95) * t83 + (Ifges(4,1) * t80 + 0.2e1 * Ifges(4,4) * t83 + t29 * t130 + (t18 + t19 + t20 + t115) * t82 + (t105 * t83 + t15 + t16 - t17) * t79) * t80; -t26 * mrSges(4,2) + (t41 - t42 + t101) * t24 + m(6) * (t24 * t38 + t93) + m(7) * (-t24 * t30 + t40 * t7 + t45 * t9) + m(5) * (-pkin(3) * t24 + t93) + t131 * (mrSges(5,3) + t107); -pkin(3) * t29 + t10 * t42 + t38 * t27 + t30 * t28 + t45 * t32 + t40 * t34 + m(7) * (t40 * t1 + t30 * t10 + t45 * t3) + (-pkin(8) * mrSges(4,2) - t67 / 0.2e1 - t66 / 0.2e1 - t65 / 0.2e1 + Ifges(4,6)) * t83 + (-t15 / 0.2e1 - t16 / 0.2e1 + t17 / 0.2e1 + t14 * mrSges(5,3) + t11 * mrSges(6,2) - t3 * mrSges(7,3) + (-Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1) * t83) * t82 + (t115 / 0.2e1 + t18 / 0.2e1 + t19 / 0.2e1 + t20 / 0.2e1 - t13 * mrSges(5,3) + t12 * mrSges(6,2) - t1 * mrSges(7,3)) * t79 + (t103 * t82 + t102 * t79 + m(5) * (-t13 * t79 + t14 * t82) + m(6) * (t11 * t82 + t12 * t79)) * pkin(9) + (Ifges(4,5) + (t49 / 0.2e1 + t50 / 0.2e1 + t51 / 0.2e1) * t82 + (t46 / 0.2e1 + t47 / 0.2e1 - t48 / 0.2e1) * t79 + (-m(5) * pkin(3) + t101) * pkin(8)) * t80 + (m(6) * t38 + t41) * t21; -0.2e1 * pkin(3) * t43 + 0.2e1 * t30 * t42 + 0.2e1 * t38 * t41 + Ifges(4,3) + m(6) * (t38 ^ 2 + t100) + m(7) * (t30 ^ 2 + t40 ^ 2 + t45 ^ 2) + m(5) * (pkin(3) ^ 2 + t100) + (t127 * t45 - t46 - t47 + t48) * t82 + (t127 * t40 + t49 + t50 + t51) * t79 + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * pkin(9) * t132; (-mrSges(5,1) + t108) * t7 + m(6) * (-pkin(4) * t7 + t4) + m(7) * (-t126 * t7 + t4) + (-mrSges(5,2) + t106) * t9; t13 * mrSges(5,1) - t12 * mrSges(6,1) - t1 * mrSges(7,1) - t14 * mrSges(5,2) + t3 * mrSges(7,2) + t11 * mrSges(6,3) - pkin(4) * t36 - t126 * t34 + (t32 + t37) * qJ(5) + m(6) * (-pkin(4) * t12 + qJ(5) * t11) + m(7) * (qJ(5) * t3 - t1 * t126) - t96 * t83 + (-Ifges(7,5) * t82 - t105 * t79) * t80 - t95; m(7) * (qJ(5) * t45 - t126 * t40) + t45 * mrSges(7,2) - t40 * mrSges(7,1) + t67 + t66 + t65 + (-pkin(4) * mrSges(6,2) + mrSges(7,3) * t126 - Ifges(7,5)) * t79 + (qJ(5) * t107 - Ifges(6,6) + Ifges(7,6)) * t82 + (m(6) * t89 - t90 - t91) * pkin(9); 0.2e1 * pkin(4) * mrSges(6,1) + 0.2e1 * t126 * mrSges(7,1) + 0.2e1 * t106 * qJ(5) + m(6) * (pkin(4) ^ 2 + t86) + m(7) * (t126 ^ 2 + t86) + t96; 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t7; m(6) * t12 + m(7) * t1 + t34 + t36; m(7) * t40 + (m(6) * pkin(9) + t107) * t79; -m(6) * pkin(4) - m(7) * t126 + t108; t133; -m(7) * t24; m(7) * t10 + t28; m(7) * t30 + t42; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
