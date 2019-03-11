% Calculate Gravitation load on the joints for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR5_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:08:12
% EndTime: 2019-03-09 07:08:14
% DurationCPUTime: 0.72s
% Computational Cost: add. (541->111), mult. (467->122), div. (0->0), fcn. (407->12), ass. (0->68)
t46 = qJ(5) + qJ(6);
t41 = sin(t46);
t50 = sin(qJ(5));
t135 = -mrSges(6,2) * t50 - mrSges(7,2) * t41;
t134 = -mrSges(6,3) - mrSges(7,3);
t42 = cos(t46);
t52 = cos(qJ(5));
t133 = t52 * mrSges(6,1) + t42 * mrSges(7,1);
t45 = pkin(11) + qJ(3);
t40 = qJ(4) + t45;
t34 = sin(t40);
t35 = cos(t40);
t132 = t35 * (-m(6) * pkin(9) + t134) + t135 * t34;
t110 = m(7) * pkin(5);
t128 = t133 * t34;
t51 = sin(qJ(1));
t53 = cos(qJ(1));
t127 = g(1) * t53 + g(2) * t51;
t126 = -t35 * mrSges(5,1) + (mrSges(5,2) + t134) * t34;
t121 = t35 * pkin(4) + t34 * pkin(9);
t37 = t52 * pkin(5) + pkin(4);
t54 = -pkin(10) - pkin(9);
t123 = -t34 * t54 + t35 * t37;
t125 = -m(6) * t121 - m(7) * t123;
t124 = t50 * t110;
t106 = pkin(4) * t34;
t38 = sin(t45);
t107 = pkin(3) * t38;
t63 = -t34 * t37 - t35 * t54;
t120 = -m(7) * (t63 - t107) - m(6) * (-t106 - t107) + t128;
t118 = mrSges(6,1) + t110;
t117 = m(5) + m(6) + m(7);
t116 = t126 + (-t133 - t135) * t35;
t49 = -pkin(7) - qJ(2);
t115 = -m(3) * qJ(2) + m(4) * t49 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t114 = t132 * t51;
t113 = t132 * t53;
t112 = m(6) * t106 - m(7) * t63 + t128;
t48 = cos(pkin(11));
t36 = t48 * pkin(2) + pkin(1);
t39 = cos(t45);
t69 = t39 * mrSges(4,1) - t38 * mrSges(4,2);
t111 = m(4) * t36 + mrSges(2,1) + m(3) * pkin(1) + t48 * mrSges(3,1) - sin(pkin(11)) * mrSges(3,2) + t69 - t126;
t88 = t53 * t42;
t94 = t51 * t41;
t5 = t35 * t94 + t88;
t89 = t53 * t41;
t93 = t51 * t42;
t6 = -t35 * t93 + t89;
t109 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t7 = -t35 * t89 + t93;
t8 = t35 * t88 + t94;
t108 = t7 * mrSges(7,1) - t8 * mrSges(7,2);
t33 = pkin(3) * t39;
t102 = g(3) * t34;
t92 = t51 * t50;
t91 = t51 * t52;
t87 = t53 * t50;
t86 = t53 * t52;
t66 = mrSges(5,1) * t34 + mrSges(5,2) * t35;
t65 = -mrSges(7,1) * t41 - mrSges(7,2) * t42;
t11 = -t35 * t87 + t91;
t9 = t35 * t92 + t86;
t44 = -pkin(8) + t49;
t18 = t33 + t36;
t12 = t35 * t86 + t92;
t10 = -t35 * t91 + t87;
t1 = [(-t92 * t110 - t12 * mrSges(6,1) - t8 * mrSges(7,1) - t11 * mrSges(6,2) - t7 * mrSges(7,2) - t117 * (t53 * t18 - t51 * t44) + t115 * t51 + (-t111 + t125) * t53) * g(2) + (-t10 * mrSges(6,1) - t6 * mrSges(7,1) - t9 * mrSges(6,2) - t5 * mrSges(7,2) + (t117 * t44 + t115 - t124) * t53 + (m(5) * t18 - m(6) * (-t18 - t121) - m(7) * (-t18 - t123) + t111) * t51) * g(1) (-g(1) * t51 + g(2) * t53) * (m(3) + m(4) + t117) (t120 * t51 + t114) * g(2) + (t120 * t53 + t113) * g(1) + (-t69 - m(5) * t33 - m(6) * (t33 + t121) - m(7) * (t33 + t123) + t116) * g(3) + t127 * (m(5) * t107 + mrSges(4,1) * t38 + mrSges(4,2) * t39 + t66) t127 * t66 + (t112 * t51 + t114) * g(2) + (t112 * t53 + t113) * g(1) + (t116 + t125) * g(3) (mrSges(6,1) * t50 + mrSges(6,2) * t52 + t124 - t65) * t102 + (-t10 * mrSges(6,2) + t118 * t9 - t109) * g(2) + (t12 * mrSges(6,2) - t11 * t118 - t108) * g(1), -g(1) * t108 - g(2) * t109 - t65 * t102];
taug  = t1(:);
