% Calculate Gravitation load on the joints for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:05:24
% EndTime: 2019-03-08 23:05:26
% DurationCPUTime: 1.06s
% Computational Cost: add. (705->113), mult. (1116->152), div. (0->0), fcn. (1278->14), ass. (0->58)
t56 = pkin(12) + qJ(6);
t52 = sin(t56);
t53 = cos(t56);
t58 = sin(pkin(12));
t60 = cos(pkin(12));
t145 = t60 * mrSges(6,1) + t53 * mrSges(7,1) - t58 * mrSges(6,2) - t52 * mrSges(7,2) + mrSges(5,1);
t141 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t50 = pkin(5) * t60 + pkin(4);
t57 = qJ(3) + qJ(4);
t54 = sin(t57);
t55 = cos(t57);
t61 = -pkin(10) - qJ(5);
t62 = sin(qJ(3));
t64 = cos(qJ(3));
t122 = m(4) * pkin(2) + t64 * mrSges(4,1) - t62 * mrSges(4,2) + mrSges(3,1) + (m(6) * pkin(4) + m(7) * t50 + t145) * t55 + (m(6) * qJ(5) - m(7) * t61 - t141) * t54;
t59 = sin(pkin(6));
t63 = sin(qJ(2));
t107 = t59 * t63;
t97 = cos(pkin(6));
t139 = -t62 * t107 + t97 * t64;
t65 = cos(qJ(2));
t95 = sin(pkin(11));
t76 = t97 * t95;
t96 = cos(pkin(11));
t41 = -t63 * t76 + t96 * t65;
t89 = t59 * t95;
t137 = -t41 * t62 + t64 * t89;
t77 = t97 * t96;
t39 = t63 * t77 + t95 * t65;
t90 = t59 * t96;
t19 = t39 * t54 + t55 * t90;
t20 = t39 * t55 - t54 * t90;
t136 = t141 * t20 + t145 * t19;
t21 = t41 * t54 - t55 * t89;
t22 = t41 * t55 + t54 * t89;
t135 = t141 * t22 + t145 * t21;
t34 = t54 * t107 - t97 * t55;
t35 = t55 * t107 + t97 * t54;
t134 = t141 * t35 + t145 * t34;
t121 = -m(4) * pkin(8) - t52 * mrSges(7,1) - t60 * mrSges(6,2) - t53 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) + (-m(7) * pkin(5) - mrSges(6,1)) * t58;
t133 = m(6) + m(7);
t126 = m(5) + t133;
t120 = -t19 * t50 - t20 * t61;
t119 = -t21 * t50 - t22 * t61;
t106 = t59 * t65;
t103 = -t34 * t50 - t35 * t61;
t87 = -t19 * pkin(4) + t20 * qJ(5);
t86 = -t21 * pkin(4) + qJ(5) * t22;
t85 = -t34 * pkin(4) + qJ(5) * t35;
t84 = t137 * pkin(3);
t78 = t139 * pkin(3);
t71 = -t39 * t62 - t64 * t90;
t69 = t71 * pkin(3);
t66 = -pkin(9) - pkin(8);
t51 = pkin(3) * t64 + pkin(2);
t40 = t96 * t63 + t65 * t76;
t38 = t95 * t63 - t65 * t77;
t1 = [(-m(2) - m(3) - m(4) - t126) * g(3) (-t126 * (-t38 * t51 - t39 * t66) + t121 * t39 + t122 * t38) * g(2) + (-t126 * (-t40 * t51 - t41 * t66) + t121 * t41 + t122 * t40) * g(1) + (-t126 * t51 * t106 + (-t122 * t65 + (t126 * t66 + t121) * t63) * t59) * g(3) (-t139 * mrSges(4,1) - (-t64 * t107 - t97 * t62) * mrSges(4,2) - m(5) * t78 - m(6) * (t78 + t85) - m(7) * (t78 + t103) + t134) * g(3) + (-t71 * mrSges(4,1) - (-t39 * t64 + t62 * t90) * mrSges(4,2) - m(5) * t69 - m(6) * (t69 + t87) - m(7) * (t69 + t120) + t136) * g(2) + (-t137 * mrSges(4,1) - (-t41 * t64 - t62 * t89) * mrSges(4,2) - m(5) * t84 - m(6) * (t84 + t86) - m(7) * (t84 + t119) + t135) * g(1) (-m(6) * t85 - m(7) * t103 + t134) * g(3) + (-m(6) * t87 - m(7) * t120 + t136) * g(2) + (-m(6) * t86 - m(7) * t119 + t135) * g(1), t133 * (-g(1) * t21 - g(2) * t19 - g(3) * t34) -g(1) * ((-t22 * t52 + t40 * t53) * mrSges(7,1) + (-t22 * t53 - t40 * t52) * mrSges(7,2)) - g(2) * ((-t20 * t52 + t38 * t53) * mrSges(7,1) + (-t20 * t53 - t38 * t52) * mrSges(7,2)) - g(3) * ((-t53 * t106 - t35 * t52) * mrSges(7,1) + (t52 * t106 - t35 * t53) * mrSges(7,2))];
taug  = t1(:);
