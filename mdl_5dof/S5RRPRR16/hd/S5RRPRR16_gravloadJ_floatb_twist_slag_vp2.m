% Calculate Gravitation load on the joints for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR16_gravloadJ_floatb_twist_slag_vp2(qJ, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:43
% EndTime: 2019-12-31 20:44:46
% DurationCPUTime: 0.78s
% Computational Cost: add. (332->85), mult. (797->119), div. (0->0), fcn. (910->10), ass. (0->50)
t53 = -m(6) * pkin(9) + mrSges(5,2) - mrSges(6,3);
t71 = m(5) + m(6);
t77 = m(4) + t71;
t85 = t77 * qJ(3);
t84 = t71 * pkin(8);
t34 = sin(qJ(5));
t38 = cos(qJ(5));
t78 = m(6) * pkin(4) + t38 * mrSges(6,1) - t34 * mrSges(6,2) + mrSges(5,1);
t35 = sin(qJ(4));
t39 = cos(qJ(4));
t66 = mrSges(4,3) - mrSges(3,2);
t83 = -t78 * t35 - t53 * t39 - t66;
t49 = t34 * mrSges(6,1) + t38 * mrSges(6,2);
t61 = mrSges(4,2) - mrSges(3,1) - mrSges(5,3);
t81 = -t49 + t61;
t80 = pkin(2) * t77 - t81 + t84;
t76 = t61 - t84;
t75 = t66 + t85;
t72 = t83 - t85;
t33 = sin(pkin(5));
t36 = sin(qJ(2));
t70 = t33 * t36;
t37 = sin(qJ(1));
t69 = t33 * t37;
t40 = cos(qJ(2));
t68 = t33 * t40;
t41 = cos(qJ(1));
t67 = t33 * t41;
t65 = pkin(2) * t68 + qJ(3) * t70;
t64 = t41 * pkin(1) + pkin(7) * t69;
t63 = cos(pkin(5));
t56 = t37 * t63;
t21 = -t36 * t56 + t40 * t41;
t60 = t21 * pkin(2) + t64;
t57 = -pkin(1) * t37 + pkin(7) * t67;
t55 = t41 * t63;
t54 = pkin(3) * t69 + t60;
t19 = t36 * t55 + t37 * t40;
t52 = -t19 * pkin(2) + t57;
t51 = mrSges(2,2) + (-mrSges(4,1) - mrSges(3,3)) * t33;
t18 = t36 * t37 - t40 * t55;
t7 = -t18 * t35 + t39 * t67;
t5 = t18 * t39 + t35 * t67;
t20 = t41 * t36 + t40 * t56;
t17 = -t35 * t68 + t63 * t39;
t4 = t20 * t35 + t39 * t69;
t3 = -t20 * t39 + t35 * t69;
t2 = t21 * t34 + t38 * t4;
t1 = t21 * t38 - t34 * t4;
t6 = [(-t41 * mrSges(2,1) - m(3) * t64 - m(4) * t60 - m(5) * t54 - t4 * mrSges(5,1) - m(6) * (pkin(4) * t4 + t54) - t2 * mrSges(6,1) - t1 * mrSges(6,2) + t53 * t3 - t75 * t20 + t51 * t37 + t76 * t21) * g(2) + (t37 * mrSges(2,1) - m(3) * t57 - m(4) * t52 + t53 * t5 + t75 * t18 - t78 * t7 + t51 * t41 + (t49 - t76) * t19 + t71 * (-pkin(3) * t67 - t52)) * g(1), (-m(4) * t65 - t71 * (pkin(8) * t68 + t65) + (t83 * t36 + t81 * t40) * t33) * g(3) + (t80 * t18 + t72 * t19) * g(2) + (t80 * t20 + t72 * t21) * g(1), t77 * (-g(1) * t20 - g(2) * t18 + g(3) * t68), (t53 * t17 - t78 * (-t63 * t35 - t39 * t68)) * g(3) + (-t78 * t5 - t53 * t7) * g(2) + (t78 * t3 + t53 * t4) * g(1), -g(1) * (mrSges(6,1) * t1 - mrSges(6,2) * t2) - g(2) * ((t19 * t38 + t34 * t7) * mrSges(6,1) + (-t19 * t34 + t38 * t7) * mrSges(6,2)) - g(3) * ((-t17 * t34 + t38 * t70) * mrSges(6,1) + (-t17 * t38 - t34 * t70) * mrSges(6,2))];
taug = t6(:);
