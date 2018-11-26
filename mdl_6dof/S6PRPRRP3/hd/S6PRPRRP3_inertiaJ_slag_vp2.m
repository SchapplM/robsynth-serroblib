% Calculate joint inertia matrix for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2018-11-23 15:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:00:42
% EndTime: 2018-11-23 15:00:43
% DurationCPUTime: 0.92s
% Computational Cost: add. (954->208), mult. (2034->289), div. (0->0), fcn. (2162->10), ass. (0->93)
t126 = Ifges(6,5) + Ifges(7,5);
t125 = Ifges(6,6) + Ifges(7,6);
t124 = Ifges(6,3) + Ifges(7,3);
t123 = -2 * mrSges(7,3);
t87 = m(6) / 0.2e1 + m(7) / 0.2e1;
t122 = 0.2e1 * t87;
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t89 = t74 ^ 2 + t77 ^ 2;
t121 = 0.2e1 * t89;
t110 = cos(qJ(4));
t70 = sin(pkin(11));
t72 = cos(pkin(11));
t75 = sin(qJ(4));
t42 = -t110 * t72 + t75 * t70;
t43 = t110 * t70 + t75 * t72;
t27 = t42 * mrSges(5,1) + t43 * mrSges(5,2);
t44 = -t72 * mrSges(4,1) + t70 * mrSges(4,2);
t120 = -t27 - t44;
t119 = t125 * t77 + t126 * t74;
t100 = t74 * mrSges(6,2);
t48 = -t77 * mrSges(6,1) + t100;
t118 = -m(6) * pkin(4) - mrSges(5,1) + t48;
t60 = t74 * mrSges(7,2);
t47 = -t77 * mrSges(7,1) + t60;
t58 = -t77 * pkin(5) - pkin(4);
t117 = m(7) * t58 + t47;
t71 = sin(pkin(6));
t76 = sin(qJ(2));
t102 = t71 * t76;
t73 = cos(pkin(6));
t36 = -t70 * t102 + t73 * t72;
t37 = t72 * t102 + t73 * t70;
t14 = -t110 * t36 + t75 * t37;
t116 = t14 ^ 2;
t97 = pkin(8) + qJ(3);
t45 = t97 * t72;
t86 = t97 * t70;
t28 = t110 * t86 + t75 * t45;
t115 = t28 ^ 2;
t67 = t72 ^ 2;
t114 = 0.2e1 * t28;
t113 = m(7) * pkin(5);
t111 = -m(4) - m(5);
t109 = Ifges(6,4) * t74;
t108 = Ifges(6,4) * t77;
t107 = Ifges(7,4) * t74;
t106 = Ifges(7,4) * t77;
t105 = t28 * t14;
t104 = t43 * t74;
t103 = t43 * t77;
t78 = cos(qJ(2));
t101 = t71 * t78;
t99 = t74 * mrSges(7,3);
t96 = -qJ(6) - pkin(9);
t22 = -t42 * mrSges(7,2) - t43 * t99;
t23 = -t42 * mrSges(6,2) - mrSges(6,3) * t104;
t95 = t22 + t23;
t24 = t42 * mrSges(7,1) - mrSges(7,3) * t103;
t25 = t42 * mrSges(6,1) - mrSges(6,3) * t103;
t94 = t24 + t25;
t57 = -t72 * pkin(3) - pkin(2);
t21 = t42 * pkin(4) - t43 * pkin(9) + t57;
t30 = t110 * t45 - t75 * t86;
t6 = t74 * t21 + t77 * t30;
t19 = mrSges(7,1) * t104 + mrSges(7,2) * t103;
t90 = t70 ^ 2 + t67;
t88 = qJ(6) * t43;
t5 = t77 * t21 - t74 * t30;
t85 = t126 * t103 + t124 * t42;
t84 = mrSges(6,1) + mrSges(7,1) + t113;
t82 = mrSges(6,1) * t74 + mrSges(6,2) * t77;
t81 = -t36 * t70 + t37 * t72;
t66 = t71 ^ 2;
t55 = t66 * t78 ^ 2;
t53 = Ifges(6,1) * t74 + t108;
t52 = Ifges(7,1) * t74 + t106;
t51 = Ifges(6,2) * t77 + t109;
t50 = Ifges(7,2) * t77 + t107;
t49 = t96 * t77;
t46 = t96 * t74;
t20 = t82 * t43;
t16 = t110 * t37 + t75 * t36;
t13 = pkin(5) * t104 + t28;
t12 = Ifges(6,5) * t42 + (Ifges(6,1) * t77 - t109) * t43;
t11 = Ifges(7,5) * t42 + (Ifges(7,1) * t77 - t107) * t43;
t10 = Ifges(6,6) * t42 + (-Ifges(6,2) * t74 + t108) * t43;
t9 = Ifges(7,6) * t42 + (-Ifges(7,2) * t74 + t106) * t43;
t8 = -t74 * t101 + t77 * t16;
t7 = -t77 * t101 - t74 * t16;
t4 = -t74 * t88 + t6;
t3 = t42 * pkin(5) - t77 * t88 + t5;
t1 = [m(2) + m(5) * (t16 ^ 2 + t116 + t55) + m(4) * (t36 ^ 2 + t37 ^ 2 + t55) + m(3) * (t66 * t76 ^ 2 + t73 ^ 2 + t55) + (t7 ^ 2 + t8 ^ 2 + t116) * t122; -t16 * t42 * mrSges(5,3) + t95 * t8 + t94 * t7 + t81 * mrSges(4,3) + (t43 * mrSges(5,3) + t19 + t20) * t14 + (-t76 * mrSges(3,2) + (mrSges(3,1) + t120) * t78) * t71 + m(6) * (t5 * t7 + t6 * t8 + t105) + m(7) * (t13 * t14 + t3 * t7 + t4 * t8) + m(5) * (-t57 * t101 + t30 * t16 + t105) + m(4) * (pkin(2) * t101 + t81 * qJ(3)); Ifges(4,2) * t67 - 0.2e1 * pkin(2) * t44 + 0.2e1 * t13 * t19 + t20 * t114 + 0.2e1 * t4 * t22 + 0.2e1 * t6 * t23 + 0.2e1 * t3 * t24 + 0.2e1 * t5 * t25 + 0.2e1 * t57 * t27 + Ifges(3,3) + (Ifges(4,1) * t70 + 0.2e1 * Ifges(4,4) * t72) * t70 + 0.2e1 * t90 * qJ(3) * mrSges(4,3) + (-0.2e1 * t30 * mrSges(5,3) + Ifges(5,2) * t42 + t85) * t42 + m(4) * (t90 * qJ(3) ^ 2 + pkin(2) ^ 2) + m(5) * (t30 ^ 2 + t57 ^ 2 + t115) + m(6) * (t5 ^ 2 + t6 ^ 2 + t115) + m(7) * (t13 ^ 2 + t3 ^ 2 + t4 ^ 2) + (mrSges(5,3) * t114 + Ifges(5,1) * t43 - 0.2e1 * Ifges(5,4) * t42 + (t11 + t12) * t77 + (-t125 * t42 - t10 - t9) * t74) * t43; t111 * t101 + (t77 * t7 + t74 * t8) * t122; -m(4) * pkin(2) + t94 * t77 + t95 * t74 + m(7) * (t77 * t3 + t74 * t4) + m(6) * (t77 * t5 + t74 * t6) + m(5) * t57 - t120; t121 * t87 - t111; -t16 * mrSges(5,2) + m(7) * (t46 * t7 - t49 * t8) + (t117 + t118) * t14 + (m(6) * pkin(9) + mrSges(6,3) + mrSges(7,3)) * (-t7 * t74 + t8 * t77); t13 * t47 - t49 * t22 + t58 * t19 + t46 * t24 - t30 * mrSges(5,2) + m(7) * (t58 * t13 + t46 * t3 - t49 * t4) - pkin(4) * t20 + t118 * t28 + (t9 / 0.2e1 + t10 / 0.2e1 + t4 * mrSges(7,3) + t6 * mrSges(6,3) + (m(6) * t6 + t23) * pkin(9)) * t77 + (t11 / 0.2e1 + t12 / 0.2e1 - t3 * mrSges(7,3) - t5 * mrSges(6,3) + (-m(6) * t5 - t25) * pkin(9)) * t74 + (Ifges(5,5) + (t52 / 0.2e1 + t53 / 0.2e1) * t77 + (-t50 / 0.2e1 - t51 / 0.2e1) * t74) * t43 + (-Ifges(5,6) + t119 / 0.2e1) * t42; m(7) * (t46 * t77 - t49 * t74); -0.2e1 * pkin(4) * t48 + 0.2e1 * t58 * t47 + Ifges(5,3) + pkin(9) * mrSges(6,3) * t121 + m(7) * (t46 ^ 2 + t49 ^ 2 + t58 ^ 2) + m(6) * (t89 * pkin(9) ^ 2 + pkin(4) ^ 2) + (t123 * t49 + t50 + t51) * t77 + (t123 * t46 + t52 + t53) * t74; (-mrSges(6,2) - mrSges(7,2)) * t8 + t84 * t7; t5 * mrSges(6,1) + t3 * mrSges(7,1) - t6 * mrSges(6,2) - t4 * mrSges(7,2) - t125 * t104 + (m(7) * t3 + t24) * pkin(5) + t85; t84 * t77 - t100 - t60; t46 * mrSges(7,1) + t49 * mrSges(7,2) - t82 * pkin(9) + (m(7) * t46 - t99) * pkin(5) + t119; (0.2e1 * mrSges(7,1) + t113) * pkin(5) + t124; m(7) * t14; m(7) * t13 + t19; 0; t117; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
