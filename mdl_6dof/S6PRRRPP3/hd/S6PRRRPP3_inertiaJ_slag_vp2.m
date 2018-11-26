% Calculate joint inertia matrix for
% S6PRRRPP3
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
% Datum: 2018-11-23 15:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:21:48
% EndTime: 2018-11-23 15:21:50
% DurationCPUTime: 1.06s
% Computational Cost: add. (696->272), mult. (1530->345), div. (0->0), fcn. (1292->8), ass. (0->104)
t133 = m(6) + m(7);
t132 = -Ifges(7,4) - Ifges(6,5);
t79 = sin(qJ(4));
t82 = cos(qJ(4));
t130 = t79 ^ 2 + t82 ^ 2;
t129 = 2 * pkin(8);
t76 = sin(pkin(6));
t84 = cos(qJ(2));
t108 = t76 * t84;
t81 = sin(qJ(2));
t109 = t76 * t81;
t77 = cos(pkin(6));
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t24 = t109 * t83 + t77 * t80;
t5 = t108 * t82 + t24 * t79;
t7 = -t108 * t79 + t24 * t82;
t127 = t5 * t79 + t7 * t82;
t22 = t109 * t80 - t77 * t83;
t21 = t22 ^ 2;
t126 = 2 * mrSges(7,1);
t125 = pkin(5) + pkin(9);
t124 = pkin(8) * t83;
t120 = Ifges(5,4) * t79;
t119 = Ifges(5,4) * t82;
t118 = Ifges(6,4) * t83;
t117 = Ifges(5,6) * t83;
t116 = Ifges(6,6) * t79;
t115 = Ifges(6,6) * t82;
t114 = Ifges(7,6) * t79;
t113 = Ifges(7,6) * t82;
t112 = qJ(5) * t7;
t111 = t22 * t80;
t110 = t24 * t83;
t107 = t79 * t80;
t106 = t80 * t82;
t105 = t82 * mrSges(7,1);
t104 = mrSges(6,1) + mrSges(7,1);
t103 = mrSges(6,2) - mrSges(7,3);
t102 = mrSges(6,3) + mrSges(7,2);
t78 = -pkin(4) - qJ(6);
t31 = -mrSges(5,1) * t83 - mrSges(5,3) * t106;
t35 = mrSges(6,1) * t106 - t83 * mrSges(6,2);
t101 = -t31 + t35;
t33 = mrSges(6,1) * t107 + mrSges(6,3) * t83;
t34 = -mrSges(7,1) * t107 - t83 * mrSges(7,2);
t100 = -t33 + t34;
t39 = -mrSges(5,1) * t82 + mrSges(5,2) * t79;
t99 = t39 - mrSges(4,1);
t98 = pkin(4) * t107 + t80 * pkin(8);
t37 = -pkin(3) * t83 - pkin(9) * t80 - pkin(2);
t13 = t82 * t124 + t79 * t37;
t32 = t83 * mrSges(7,3) + t80 * t105;
t97 = t130 * pkin(9) ^ 2;
t96 = qJ(5) * t82;
t94 = Ifges(6,1) + Ifges(5,3) + Ifges(7,1);
t93 = -qJ(5) * t79 - pkin(3);
t58 = t79 * t124;
t12 = t37 * t82 - t58;
t92 = t132 * t107 + (-Ifges(5,5) - Ifges(7,5)) * t106;
t10 = qJ(5) * t83 - t13;
t90 = t79 * mrSges(5,1) + t82 * mrSges(5,2);
t89 = -t79 * mrSges(6,2) - t82 * mrSges(6,3);
t88 = t127 * pkin(9);
t87 = pkin(8) ^ 2;
t85 = qJ(5) ^ 2;
t75 = t83 ^ 2;
t73 = t80 ^ 2;
t71 = t76 ^ 2;
t70 = t83 * pkin(4);
t67 = t73 * t87;
t65 = Ifges(5,5) * t79;
t64 = Ifges(7,5) * t79;
t63 = Ifges(5,6) * t82;
t60 = t71 * t84 ^ 2;
t49 = t125 * t82;
t48 = t125 * t79;
t47 = Ifges(5,1) * t79 + t119;
t46 = Ifges(5,2) * t82 + t120;
t45 = -Ifges(6,2) * t79 - t115;
t44 = -Ifges(7,2) * t82 + t114;
t43 = -Ifges(6,3) * t82 - t116;
t42 = Ifges(7,3) * t79 - t113;
t41 = -mrSges(7,2) * t79 - mrSges(7,3) * t82;
t40 = -mrSges(4,1) * t83 + mrSges(4,2) * t80;
t38 = mrSges(6,2) * t82 - mrSges(6,3) * t79;
t36 = -pkin(4) * t82 + t93;
t30 = mrSges(5,2) * t83 - mrSges(5,3) * t107;
t28 = t90 * t80;
t27 = t89 * t80;
t26 = (-mrSges(7,2) * t82 + mrSges(7,3) * t79) * t80;
t25 = t78 * t82 + t93;
t20 = -t80 * t96 + t98;
t19 = -t118 + (-Ifges(6,2) * t82 + t116) * t80;
t18 = -Ifges(7,4) * t83 + (Ifges(7,2) * t79 + t113) * t80;
t17 = -Ifges(6,5) * t83 + (Ifges(6,3) * t79 - t115) * t80;
t16 = -Ifges(7,5) * t83 + (Ifges(7,3) * t82 + t114) * t80;
t15 = -Ifges(5,5) * t83 + (Ifges(5,1) * t82 - t120) * t80;
t14 = -t117 + (-Ifges(5,2) * t79 + t119) * t80;
t11 = -t12 + t70;
t9 = (qJ(6) * t79 - t96) * t80 + t98;
t2 = -pkin(5) * t107 - t10;
t1 = qJ(6) * t83 + t58 + t70 + (pkin(5) * t80 - t37) * t82;
t3 = [m(2) + m(4) * (t24 ^ 2 + t21 + t60) + m(3) * (t71 * t81 ^ 2 + t77 ^ 2 + t60) + (m(5) + t133) * (t5 ^ 2 + t7 ^ 2 + t21); mrSges(4,3) * t110 + (-t81 * mrSges(3,2) + (mrSges(3,1) - t40) * t84) * t76 + (t30 + t100) * t7 + (t32 + t101) * t5 + (t80 * mrSges(4,3) + t26 + t27 + t28) * t22 + m(5) * (pkin(8) * t111 - t12 * t5 + t13 * t7) + m(6) * (-t10 * t7 + t11 * t5 + t20 * t22) + m(7) * (t1 * t5 + t2 * t7 + t22 * t9) + m(4) * (pkin(2) * t108 + (t110 + t111) * pkin(8)); -0.2e1 * pkin(2) * t40 + 0.2e1 * t1 * t32 + 0.2e1 * t10 * t33 + 0.2e1 * t11 * t35 + 0.2e1 * t12 * t31 + 0.2e1 * t13 * t30 + 0.2e1 * t2 * t34 + 0.2e1 * t20 * t27 + 0.2e1 * t9 * t26 + Ifges(3,3) + (t73 + t75) * mrSges(4,3) * t129 + m(4) * (pkin(2) ^ 2 + t75 * t87 + t67) + m(5) * (t12 ^ 2 + t13 ^ 2 + t67) + m(6) * (t10 ^ 2 + t11 ^ 2 + t20 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) + ((Ifges(4,2) + t94) * t83 + t92) * t83 + (Ifges(4,1) * t80 + 0.2e1 * Ifges(4,4) * t83 + t28 * t129 + (t15 + t16 - t19 + t118) * t82 + (-t14 + t17 + t18 + t117) * t79) * t80; -t24 * mrSges(4,2) + (t38 + t41 + t99) * t22 + m(5) * (-pkin(3) * t22 + t88) + m(6) * (t22 * t36 + t88) + m(7) * (t22 * t25 + t48 * t5 + t49 * t7) + t127 * (mrSges(5,3) + t104); -pkin(3) * t28 + t25 * t26 + t36 * t27 + t48 * t32 + t49 * t34 + t9 * t41 + m(7) * (t1 * t48 + t2 * t49 + t25 * t9) + (-pkin(8) * mrSges(4,2) - t64 / 0.2e1 - t65 / 0.2e1 - t63 / 0.2e1 + Ifges(4,6)) * t83 + (t14 / 0.2e1 - t17 / 0.2e1 - t18 / 0.2e1 + t13 * mrSges(5,3) - t10 * mrSges(6,1) + t2 * mrSges(7,1) + (Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1) * t83) * t82 + (t118 / 0.2e1 + t15 / 0.2e1 + t16 / 0.2e1 - t19 / 0.2e1 - t12 * mrSges(5,3) + t11 * mrSges(6,1) + t1 * mrSges(7,1)) * t79 + ((t30 - t33) * t82 + t101 * t79 + m(5) * (-t12 * t79 + t13 * t82) + m(6) * (-t10 * t82 + t11 * t79)) * pkin(9) + (Ifges(4,5) + (-t45 / 0.2e1 + t47 / 0.2e1 + t42 / 0.2e1) * t82 + (t43 / 0.2e1 + t44 / 0.2e1 - t46 / 0.2e1) * t79 + (-m(5) * pkin(3) + t99) * pkin(8)) * t80 + (m(6) * t36 + t38) * t20; -0.2e1 * pkin(3) * t39 + 0.2e1 * t25 * t41 + 0.2e1 * t36 * t38 + Ifges(4,3) + m(6) * (t36 ^ 2 + t97) + m(7) * (t25 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (pkin(3) ^ 2 + t97) + (t126 * t49 - t43 - t44 + t46) * t82 + (t126 * t48 + t42 - t45 + t47) * t79 + 0.2e1 * (mrSges(6,1) + mrSges(5,3)) * pkin(9) * t130; (-mrSges(5,2) + t102) * t7 + (-mrSges(5,1) + t103) * t5 + m(6) * (-pkin(4) * t5 + t112) + m(7) * (t5 * t78 + t112); t12 * mrSges(5,1) - t13 * mrSges(5,2) + t11 * mrSges(6,2) + t2 * mrSges(7,2) - t10 * mrSges(6,3) - t1 * mrSges(7,3) - pkin(4) * t35 + t78 * t32 + (-Ifges(6,4) * t82 - Ifges(5,6) * t79) * t80 + t100 * qJ(5) + m(7) * (qJ(5) * t2 + t1 * t78) + m(6) * (-pkin(4) * t11 - qJ(5) * t10) - t94 * t83 - t92; t65 + t63 + t64 + m(7) * (qJ(5) * t49 + t48 * t78) + t49 * mrSges(7,2) - t48 * mrSges(7,3) + (-pkin(4) * mrSges(6,1) + t78 * mrSges(7,1) - Ifges(6,4)) * t79 + (qJ(5) * t104 + t132) * t82 + (m(6) * (-pkin(4) * t79 + t96) - t89 - t90) * pkin(9); -0.2e1 * pkin(4) * mrSges(6,2) - 0.2e1 * t78 * mrSges(7,3) + 0.2e1 * t102 * qJ(5) + m(7) * (t78 ^ 2 + t85) + m(6) * (pkin(4) ^ 2 + t85) + t94; 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t5; m(6) * t11 + m(7) * t1 + t32 + t35; m(7) * t48 + (m(6) * pkin(9) + t104) * t79; -m(6) * pkin(4) + m(7) * t78 + t103; t133; m(7) * t7; m(7) * t2 + t34; m(7) * t49 + t105; m(7) * qJ(5) + mrSges(7,2); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
