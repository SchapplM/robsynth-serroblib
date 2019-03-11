% Calculate Gravitation load on the joints for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRRR4_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR4_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:46:45
% EndTime: 2019-03-10 03:46:47
% DurationCPUTime: 1.04s
% Computational Cost: add. (730->152), mult. (717->181), div. (0->0), fcn. (676->12), ass. (0->81)
t122 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t54 = qJ(3) + qJ(4);
t46 = cos(t54);
t40 = pkin(4) * t46;
t58 = cos(qJ(3));
t50 = t58 * pkin(3);
t36 = t40 + t50;
t47 = qJ(5) + t54;
t42 = cos(t47);
t37 = pkin(5) * t42;
t27 = t37 + t36;
t21 = pkin(2) + t27;
t34 = pkin(2) + t36;
t44 = qJ(6) + t47;
t38 = sin(t44);
t39 = cos(t44);
t41 = sin(t47);
t43 = t50 + pkin(2);
t45 = sin(t54);
t55 = sin(qJ(3));
t121 = -m(4) * pkin(2) - m(5) * t43 - m(6) * t34 - m(7) * t21 - t58 * mrSges(4,1) - t46 * mrSges(5,1) - t42 * mrSges(6,1) - t39 * mrSges(7,1) + t55 * mrSges(4,2) + t45 * mrSges(5,2) + t41 * mrSges(6,2) + t38 * mrSges(7,2);
t61 = -pkin(9) - pkin(8);
t53 = -pkin(10) + t61;
t48 = -pkin(11) + t53;
t120 = -m(4) * pkin(8) + m(5) * t61 + m(6) * t53 + m(7) * t48 - t122;
t57 = sin(qJ(1));
t60 = cos(qJ(1));
t119 = g(1) * t60 + g(2) * t57;
t106 = m(5) * pkin(3);
t73 = -mrSges(7,1) * t38 - mrSges(7,2) * t39;
t118 = mrSges(6,1) * t41 + mrSges(6,2) * t42 - t73;
t117 = mrSges(4,1) + t106;
t59 = cos(qJ(2));
t87 = t59 * t60;
t7 = -t38 * t87 + t39 * t57;
t8 = t38 * t57 + t39 * t87;
t102 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t15 = -t41 * t87 + t42 * t57;
t16 = t41 * t57 + t42 * t87;
t116 = -t15 * mrSges(6,1) + t16 * mrSges(6,2) - t102;
t88 = t57 * t59;
t5 = t38 * t88 + t39 * t60;
t6 = t38 * t60 - t39 * t88;
t103 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t13 = t41 * t88 + t42 * t60;
t14 = t41 * t60 - t42 * t88;
t115 = t13 * mrSges(6,1) - t14 * mrSges(6,2) - t103;
t114 = mrSges(5,1) * t45 + mrSges(5,2) * t46 + t118;
t24 = -t45 * t87 + t46 * t57;
t25 = t45 * t57 + t46 * t87;
t113 = -t24 * mrSges(5,1) + t25 * mrSges(5,2) + t116;
t22 = t45 * t88 + t46 * t60;
t23 = t45 * t60 - t46 * t88;
t112 = t22 * mrSges(5,1) - t23 * mrSges(5,2) + t115;
t111 = -m(3) - m(4) - m(5) - m(6) - m(7);
t101 = pkin(3) * t55;
t100 = pkin(4) * t45;
t99 = pkin(5) * t41;
t32 = -t99 - t100;
t26 = -t32 + t101;
t35 = t100 + t101;
t109 = -m(6) * t35 - m(7) * t26;
t56 = sin(qJ(2));
t77 = t59 * mrSges(3,1) - t56 * mrSges(3,2);
t108 = t122 * t56 + mrSges(2,1) + t77;
t107 = mrSges(2,2) - mrSges(3,3) + t109;
t105 = m(6) * pkin(4);
t104 = m(7) * pkin(5);
t96 = g(3) * t56;
t94 = t55 * t57;
t93 = t55 * t60;
t78 = pkin(2) * t59 + pkin(8) * t56;
t72 = t21 * t59 - t48 * t56;
t71 = t34 * t59 - t53 * t56;
t70 = t59 * t43 - t56 * t61;
t30 = -t55 * t87 + t57 * t58;
t28 = t55 * t88 + t58 * t60;
t33 = t37 + t40;
t31 = t58 * t87 + t94;
t29 = -t58 * t88 + t93;
t1 = [(-t94 * t106 - t31 * mrSges(4,1) - t25 * mrSges(5,1) - t16 * mrSges(6,1) - t8 * mrSges(7,1) - t30 * mrSges(4,2) - t24 * mrSges(5,2) - t15 * mrSges(6,2) - t7 * mrSges(7,2) + t111 * (t60 * pkin(1) + t57 * pkin(7)) + t107 * t57 + (-m(4) * t78 - m(5) * t70 - m(6) * t71 - m(7) * t72 - t108) * t60) * g(2) + (-t93 * t106 - t29 * mrSges(4,1) - t23 * mrSges(5,1) - t14 * mrSges(6,1) - t6 * mrSges(7,1) - t28 * mrSges(4,2) - t22 * mrSges(5,2) - t13 * mrSges(6,2) - t5 * mrSges(7,2) + (m(3) * pkin(1) - m(4) * (-pkin(1) - t78) - m(5) * (-pkin(1) - t70) - m(6) * (-pkin(1) - t71) - m(7) * (-pkin(1) - t72) + t108) * t57 + (t111 * pkin(7) + t107) * t60) * g(1), -g(3) * t77 + (t121 * g(3) + t119 * (mrSges(3,2) + t120)) * t59 + (t120 * g(3) + t119 * (mrSges(3,1) - t121)) * t56 (m(5) * t101 + mrSges(4,1) * t55 + mrSges(4,2) * t58 - t109 + t114) * t96 + (-t29 * mrSges(4,2) - m(6) * (-t35 * t88 - t36 * t60) - m(7) * (-t26 * t88 - t27 * t60) + t117 * t28 + t112) * g(2) + (t31 * mrSges(4,2) - m(6) * (-t35 * t87 + t36 * t57) - m(7) * (-t26 * t87 + t27 * t57) - t117 * t30 + t113) * g(1) (m(6) * t100 - m(7) * t32 + t114) * t96 + (t22 * t105 - m(7) * (t32 * t88 - t33 * t60) + t112) * g(2) + (-t24 * t105 - m(7) * (t32 * t87 + t33 * t57) + t113) * g(1) (m(7) * t99 + t118) * t96 + (t13 * t104 + t115) * g(2) + (-t15 * t104 + t116) * g(1), -g(1) * t102 - g(2) * t103 - t73 * t96];
taug  = t1(:);
