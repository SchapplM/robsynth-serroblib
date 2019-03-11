% Calculate potential energy for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:44
% EndTime: 2019-03-08 20:26:44
% DurationCPUTime: 0.78s
% Computational Cost: add. (421->91), mult. (867->97), div. (0->0), fcn. (1059->14), ass. (0->47)
t45 = cos(qJ(2));
t80 = t45 * mrSges(3,2);
t41 = sin(qJ(4));
t44 = cos(qJ(4));
t34 = qJ(5) + qJ(6);
t31 = sin(t34);
t32 = cos(t34);
t40 = sin(qJ(5));
t43 = cos(qJ(5));
t69 = -m(6) * pkin(4) - m(7) * (pkin(5) * t43 + pkin(4)) - t43 * mrSges(6,1) - t32 * mrSges(7,1) + t40 * mrSges(6,2) + t31 * mrSges(7,2) - mrSges(5,1);
t72 = -m(6) * pkin(9) + m(7) * (-pkin(10) - pkin(9)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t79 = t41 * t72 + t44 * t69 - mrSges(4,1);
t37 = sin(pkin(6));
t39 = cos(pkin(6));
t42 = sin(qJ(2));
t64 = t39 * t42;
t76 = -mrSges(3,3) - mrSges(4,3);
t78 = -t64 * mrSges(3,1) - t39 * t80 - mrSges(2,2) + (m(3) * pkin(7) - t41 * t69 + t44 * t72 - t76) * t37;
t77 = -m(5) - m(6);
t35 = sin(pkin(12));
t61 = cos(pkin(12));
t21 = -t42 * t35 + t45 * t61;
t75 = -m(3) - m(1) - m(2);
t71 = -m(3) * pkin(1) - t45 * mrSges(3,1) + t42 * mrSges(3,2) - mrSges(2,1);
t70 = m(7) * (pkin(5) * t40 + pkin(8)) + t40 * mrSges(6,1) + t31 * mrSges(7,1) + t43 * mrSges(6,2) + t32 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3);
t60 = qJ(1) + r_base(3);
t19 = pkin(2) * t64 + (-pkin(7) - qJ(3)) * t37;
t29 = pkin(2) * t45 + pkin(1);
t36 = sin(pkin(11));
t38 = cos(pkin(11));
t57 = t38 * t19 + t36 * t29 + r_base(2);
t56 = t39 * pkin(7) + t60;
t20 = -t45 * t35 - t42 * t61;
t18 = t20 * t39;
t8 = -t18 * t38 + t21 * t36;
t55 = t8 * pkin(3) + t57;
t54 = -t36 * t19 + t38 * t29 + r_base(1);
t53 = t37 * t42 * pkin(2) + t39 * qJ(3) + t56;
t10 = t18 * t36 + t21 * t38;
t52 = t10 * pkin(3) + t54;
t17 = t20 * t37;
t51 = -t17 * pkin(3) + t53;
t48 = t39 * t21;
t16 = t21 * t37;
t9 = t20 * t38 - t36 * t48;
t7 = t36 * t20 + t38 * t48;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t60 - mrSges(2,3) - m(3) * t56 - (t42 * mrSges(3,1) + t80) * t37 - m(4) * t53 + t17 * mrSges(4,1) - m(7) * t51 + t77 * (-t16 * pkin(8) + t51) + t76 * t39 + t69 * (-t17 * t44 + t39 * t41) + t70 * t16 + t72 * (-t17 * t41 - t39 * t44)) * g(3) + (-m(4) * t57 - m(7) * t55 - mrSges(1,2) + t71 * t36 + t75 * r_base(2) + t77 * (-t7 * pkin(8) + t55) + t70 * t7 + t79 * t8 + t78 * t38) * g(2) + (-m(4) * t54 - m(7) * t52 - mrSges(1,1) + t71 * t38 + t75 * r_base(1) + t77 * (-t9 * pkin(8) + t52) + t70 * t9 + t79 * t10 - t78 * t36) * g(1);
U  = t1;
