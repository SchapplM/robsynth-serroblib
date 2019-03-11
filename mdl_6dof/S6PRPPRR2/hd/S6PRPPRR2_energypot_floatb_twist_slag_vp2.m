% Calculate potential energy for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:35
% EndTime: 2019-03-08 19:17:35
% DurationCPUTime: 0.74s
% Computational Cost: add. (357->92), mult. (738->107), div. (0->0), fcn. (881->12), ass. (0->49)
t44 = cos(qJ(2));
t78 = t44 * mrSges(3,2);
t77 = -m(6) - m(7);
t76 = mrSges(4,2) - mrSges(5,3);
t33 = sin(pkin(11));
t36 = cos(pkin(11));
t41 = sin(qJ(2));
t75 = t41 * t33 - t36 * t44;
t74 = -m(3) - m(1) - m(2);
t73 = -mrSges(3,3) - mrSges(4,3) - mrSges(5,1);
t72 = m(7) * pkin(9) - mrSges(6,2) + mrSges(7,3);
t71 = m(3) * pkin(7) - t73;
t38 = cos(pkin(6));
t61 = t38 * t41;
t70 = t61 * mrSges(3,1) + t38 * t78 + mrSges(2,2);
t39 = sin(qJ(6));
t42 = cos(qJ(6));
t69 = -t39 * mrSges(7,1) - t42 * mrSges(7,2) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3);
t68 = -m(3) * pkin(1) - t44 * mrSges(3,1) + t41 * mrSges(3,2) - mrSges(2,1);
t67 = -m(7) * pkin(5) - t42 * mrSges(7,1) + t39 * mrSges(7,2) - mrSges(6,1);
t34 = sin(pkin(10));
t35 = sin(pkin(6));
t66 = t34 * t35;
t37 = cos(pkin(10));
t65 = t35 * t37;
t40 = sin(qJ(5));
t64 = t35 * t40;
t43 = cos(qJ(5));
t63 = t35 * t43;
t58 = qJ(1) + r_base(3);
t21 = pkin(2) * t61 + (-pkin(7) - qJ(3)) * t35;
t29 = pkin(2) * t44 + pkin(1);
t57 = t37 * t21 + t34 * t29 + r_base(2);
t56 = t38 * pkin(7) + t58;
t55 = t33 * t44 + t41 * t36;
t54 = -t21 * t34 + t37 * t29 + r_base(1);
t53 = t35 * t41 * pkin(2) + t38 * qJ(3) + t56;
t52 = t55 * t38;
t51 = t75 * t38;
t8 = -t34 * t55 - t37 * t51;
t9 = -t34 * t75 + t37 * t52;
t50 = t9 * pkin(3) - qJ(4) * t8 + t57;
t10 = t34 * t51 - t37 * t55;
t11 = -t34 * t52 - t37 * t75;
t49 = t11 * pkin(3) - qJ(4) * t10 + t54;
t19 = t75 * t35;
t20 = t55 * t35;
t48 = t20 * pkin(3) + qJ(4) * t19 + t53;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t58 - mrSges(2,3) - m(3) * t56 - (t41 * mrSges(3,1) + t78) * t35 - m(4) * t53 - m(5) * t48 + t77 * (t38 * pkin(4) + pkin(8) * t20 + t48) + t76 * t19 - t72 * (-t19 * t43 + t38 * t40) + t73 * t38 + t67 * (t19 * t40 + t38 * t43) + t69 * t20) * g(3) + (-m(4) * t57 - m(5) * t50 - mrSges(1,2) - t70 * t37 + t68 * t34 + t74 * r_base(2) - t76 * t8 + t77 * (-pkin(4) * t65 + pkin(8) * t9 + t50) + t72 * (t37 * t64 - t43 * t8) + t67 * (-t37 * t63 - t40 * t8) + t69 * t9 + t71 * t65) * g(2) + (-m(4) * t54 - m(5) * t49 - mrSges(1,1) + t70 * t34 + t68 * t37 + t74 * r_base(1) + t77 * (pkin(4) * t66 + pkin(8) * t11 + t49) - t76 * t10 - t72 * (t10 * t43 + t34 * t64) - t71 * t66 + t67 * (-t10 * t40 + t34 * t63) + t69 * t11) * g(1);
U  = t1;
