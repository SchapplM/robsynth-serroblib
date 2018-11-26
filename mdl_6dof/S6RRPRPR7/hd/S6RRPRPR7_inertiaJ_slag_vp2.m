% Calculate joint inertia matrix for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:04:34
% EndTime: 2018-11-23 17:04:35
% DurationCPUTime: 1.09s
% Computational Cost: add. (1633->239), mult. (2864->326), div. (0->0), fcn. (3020->8), ass. (0->85)
t62 = sin(pkin(10));
t63 = cos(pkin(10));
t65 = sin(qJ(4));
t68 = cos(qJ(4));
t32 = t62 * t65 - t63 * t68;
t34 = t62 * t68 + t63 * t65;
t64 = sin(qJ(6));
t67 = cos(qJ(6));
t88 = t64 ^ 2 + t67 ^ 2;
t84 = t88 * mrSges(7,3);
t43 = -t67 * mrSges(7,1) + mrSges(7,2) * t64;
t91 = mrSges(6,1) - t43;
t116 = t68 * mrSges(5,1) - t65 * mrSges(5,2) - (mrSges(6,2) - t84) * t34 - t91 * t32;
t115 = m(6) * pkin(4);
t66 = sin(qJ(2));
t69 = cos(qJ(2));
t113 = t66 ^ 2 + t69 ^ 2;
t49 = pkin(4) * t62 + pkin(9);
t82 = t88 * t49;
t35 = -t65 * t66 - t68 * t69;
t36 = -t65 * t69 + t66 * t68;
t17 = -t63 * t35 + t36 * t62;
t18 = t35 * t62 + t36 * t63;
t99 = t18 * t64;
t10 = -mrSges(7,2) * t17 - mrSges(7,3) * t99;
t98 = t18 * t67;
t11 = mrSges(7,1) * t17 - mrSges(7,3) * t98;
t112 = t67 * t10 - t64 * t11;
t111 = -m(4) * pkin(2) - mrSges(4,1);
t103 = pkin(7) - pkin(8);
t85 = t103 * t69;
t86 = t103 * t66;
t26 = t65 * t86 + t68 * t85;
t13 = qJ(5) * t35 + t26;
t25 = -t65 * t85 + t68 * t86;
t74 = -t36 * qJ(5) + t25;
t6 = t13 * t62 - t63 * t74;
t109 = t6 ^ 2;
t108 = 0.2e1 * t6;
t107 = t32 ^ 2;
t42 = -t69 * pkin(2) - t66 * qJ(3) - pkin(1);
t30 = t69 * pkin(3) - t42;
t24 = -pkin(4) * t35 + t30;
t106 = 0.2e1 * t24;
t105 = -0.2e1 * t42;
t104 = t67 / 0.2e1;
t102 = t32 * t6;
t101 = Ifges(7,4) * t64;
t100 = Ifges(7,4) * t67;
t70 = -pkin(2) - pkin(3);
t40 = -qJ(3) * t65 + t68 * t70;
t39 = -pkin(4) + t40;
t41 = qJ(3) * t68 + t65 * t70;
t22 = t39 * t63 - t41 * t62;
t97 = t22 * mrSges(6,1);
t23 = t62 * t39 + t63 * t41;
t96 = t23 * mrSges(6,2);
t95 = t40 * mrSges(5,1);
t94 = t41 * mrSges(5,2);
t90 = Ifges(7,5) * t98 + Ifges(7,3) * t17;
t89 = t113 * pkin(7) ^ 2;
t44 = Ifges(7,5) * t64 + Ifges(7,6) * t67;
t87 = t44 / 0.2e1 - Ifges(6,6);
t21 = -pkin(9) + t23;
t83 = t88 * t21;
t5 = pkin(5) * t17 - pkin(9) * t18 + t24;
t8 = t63 * t13 + t62 * t74;
t1 = t5 * t67 - t64 * t8;
t2 = t5 * t64 + t67 * t8;
t80 = -t1 * t64 + t2 * t67;
t78 = t63 * mrSges(6,1) - t62 * mrSges(6,2);
t77 = mrSges(7,1) * t64 + mrSges(7,2) * t67;
t45 = Ifges(7,2) * t67 + t101;
t46 = Ifges(7,1) * t64 + t100;
t76 = t67 * t45 + t64 * t46 + Ifges(5,3) + Ifges(6,3);
t75 = t46 * t104 - t64 * t45 / 0.2e1 + Ifges(6,5);
t73 = t25 * mrSges(5,1) - t26 * mrSges(5,2) - t8 * mrSges(6,2) + Ifges(5,5) * t36 + Ifges(5,6) * t35;
t50 = -pkin(4) * t63 - pkin(5);
t31 = t34 ^ 2;
t20 = pkin(5) - t22;
t15 = t18 * mrSges(6,2);
t9 = t77 * t18;
t4 = Ifges(7,5) * t17 + (Ifges(7,1) * t67 - t101) * t18;
t3 = Ifges(7,6) * t17 + (-Ifges(7,2) * t64 + t100) * t18;
t7 = [0.2e1 * t30 * (-t35 * mrSges(5,1) + t36 * mrSges(5,2)) + t35 * (Ifges(5,4) * t36 + Ifges(5,2) * t35) + t36 * (Ifges(5,1) * t36 + Ifges(5,4) * t35) + t15 * t106 + t9 * t108 + 0.2e1 * t2 * t10 + 0.2e1 * t1 * t11 + Ifges(2,3) + (0.2e1 * pkin(1) * mrSges(3,1) + mrSges(4,1) * t105 + (Ifges(4,3) + Ifges(3,2)) * t69) * t69 + (mrSges(6,1) * t106 - 0.2e1 * t8 * mrSges(6,3) + Ifges(6,2) * t17 + t90) * t17 + (mrSges(6,3) * t108 + Ifges(6,1) * t18 - t64 * t3 + t67 * t4 + (-Ifges(7,6) * t64 - (2 * Ifges(6,4))) * t17) * t18 + m(4) * (t42 ^ 2 + t89) + m(3) * (pkin(1) ^ 2 + t89) + m(5) * (t25 ^ 2 + t26 ^ 2 + t30 ^ 2) + m(6) * (t24 ^ 2 + t8 ^ 2 + t109) + m(7) * (t1 ^ 2 + t2 ^ 2 + t109) + 0.2e1 * (-t25 * t36 + t26 * t35) * mrSges(5,3) + 0.2e1 * (mrSges(4,2) + mrSges(3,3)) * pkin(7) * t113 + (-0.2e1 * pkin(1) * mrSges(3,2) + mrSges(4,3) * t105 + 0.2e1 * (Ifges(3,4) - Ifges(4,5)) * t69 + (Ifges(4,1) + Ifges(3,1)) * t66) * t66; t20 * t9 + t91 * t6 + (t41 * t35 - t40 * t36) * mrSges(5,3) + (qJ(3) * mrSges(4,2) + Ifges(3,6) - Ifges(4,6)) * t69 + (-t2 * mrSges(7,3) + t21 * t10 - t3 / 0.2e1) * t67 + (-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t66 + (t1 * mrSges(7,3) - t21 * t11 - t4 / 0.2e1) * t64 + (-t23 * mrSges(6,3) - t87) * t17 + m(7) * (t20 * t6 + t80 * t21) + m(5) * (t25 * t40 + t26 * t41) + m(6) * (-t22 * t6 + t23 * t8) + (-t22 * mrSges(6,3) - t75) * t18 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t69 + (-mrSges(3,1) + t111) * t66) * pkin(7) - t73; 0.2e1 * pkin(2) * mrSges(4,1) - 0.2e1 * t95 - 0.2e1 * t97 + 0.2e1 * t94 + 0.2e1 * t96 + 0.2e1 * qJ(3) * mrSges(4,3) - 0.2e1 * t20 * t43 + Ifges(4,2) + Ifges(3,3) - 0.2e1 * t21 * t84 + m(6) * (t22 ^ 2 + t23 ^ 2) + m(7) * (t88 * t21 ^ 2 + t20 ^ 2) + m(5) * (t40 ^ 2 + t41 ^ 2) + m(4) * (pkin(2) ^ 2 + qJ(3) ^ 2) + t76; (m(4) * pkin(7) + mrSges(4,2)) * t66 + (t18 * mrSges(6,3) + t9) * t32 + (t35 * t65 - t36 * t68) * mrSges(5,3) + (-t17 * mrSges(6,3) + t112) * t34 + m(7) * (t80 * t34 + t102) + m(6) * (t34 * t8 + t102) + m(5) * (t25 * t68 + t26 * t65); m(6) * (-t22 * t32 + t23 * t34) + m(7) * (t20 * t32 + t34 * t83) + m(5) * (t40 * t68 + t41 * t65) + t111 - t116; m(4) + m(5) * (t65 ^ 2 + t68 ^ 2) + m(6) * (t31 + t107) + m(7) * (t88 * t31 + t107); t64 * t4 / 0.2e1 + t3 * t104 - t6 * mrSges(6,1) + t6 * t43 + t87 * t17 + t80 * mrSges(7,3) + t75 * t18 + (m(6) * (-t6 * t63 + t62 * t8) + (-t17 * t62 - t18 * t63) * mrSges(6,3)) * pkin(4) + t73 + (m(7) * t6 + t9) * t50 + (m(7) * t80 + t112) * t49; m(7) * (t20 * t50 + t21 * t82) - t94 + t95 - t96 + t97 + (-t50 + t20) * t43 + (m(6) * (t22 * t63 + t23 * t62) - t78) * pkin(4) + (t83 - t82) * mrSges(7,3) - t76; m(7) * (t50 * t32 + t34 * t82) + (-t32 * t63 + t34 * t62) * t115 + t116; 0.2e1 * t50 * t43 + m(7) * (t88 * t49 ^ 2 + t50 ^ 2) + t76 + 0.2e1 * mrSges(7,3) * t82 + (0.2e1 * t78 + (t62 ^ 2 + t63 ^ 2) * t115) * pkin(4); t17 * mrSges(6,1) + t64 * t10 + t67 * t11 + t15 + m(7) * (t1 * t67 + t2 * t64) + m(6) * t24; 0; 0; 0; m(7) * t88 + m(6); mrSges(7,1) * t1 - mrSges(7,2) * t2 - Ifges(7,6) * t99 + t90; -t77 * t21 - t44; -t77 * t34; -t77 * t49 + t44; -t43; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t7(1) t7(2) t7(4) t7(7) t7(11) t7(16); t7(2) t7(3) t7(5) t7(8) t7(12) t7(17); t7(4) t7(5) t7(6) t7(9) t7(13) t7(18); t7(7) t7(8) t7(9) t7(10) t7(14) t7(19); t7(11) t7(12) t7(13) t7(14) t7(15) t7(20); t7(16) t7(17) t7(18) t7(19) t7(20) t7(21);];
Mq  = res;
