% Calculate potential energy for
% S6PRPRRR1
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:14
% EndTime: 2019-03-08 20:22:14
% DurationCPUTime: 0.84s
% Computational Cost: add. (404->99), mult. (738->114), div. (0->0), fcn. (880->14), ass. (0->57)
t50 = cos(qJ(4));
t86 = t50 * mrSges(5,2);
t51 = cos(qJ(2));
t85 = t51 * mrSges(3,2);
t84 = -m(4) - m(5);
t83 = -m(6) - m(7);
t82 = mrSges(3,3) + mrSges(4,3);
t81 = -m(1) - m(2) - m(3);
t80 = m(3) * pkin(7) + t82;
t79 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t48 = sin(qJ(2));
t78 = -m(3) * pkin(1) - t51 * mrSges(3,1) + t48 * mrSges(3,2) - mrSges(2,1);
t47 = sin(qJ(4));
t77 = -m(5) * pkin(3) - t50 * mrSges(5,1) + t47 * mrSges(5,2) - mrSges(4,1);
t46 = sin(qJ(6));
t49 = cos(qJ(6));
t76 = -m(7) * pkin(5) - t49 * mrSges(7,1) + t46 * mrSges(7,2) - mrSges(6,1);
t42 = sin(pkin(6));
t45 = cos(pkin(6));
t66 = t45 * t48;
t75 = t66 * mrSges(3,1) - t42 * t86 + t45 * t85 + mrSges(2,2);
t74 = m(5) * pkin(8) + t46 * mrSges(7,1) + t49 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t41 = sin(pkin(11));
t73 = t41 * t42;
t44 = cos(pkin(11));
t72 = t42 * t44;
t71 = t42 * t47;
t70 = t42 * t48;
t43 = cos(pkin(12));
t68 = t43 * t51;
t67 = t45 * t47;
t64 = t41 * t71;
t63 = t44 * t71;
t62 = qJ(1) + r_base(3);
t22 = pkin(2) * t66 + (-pkin(7) - qJ(3)) * t42;
t34 = pkin(2) * t51 + pkin(1);
t61 = t44 * t22 + t41 * t34 + r_base(2);
t60 = t45 * pkin(7) + t62;
t40 = sin(pkin(12));
t59 = t40 * t51 + t43 * t48;
t24 = -t40 * t48 + t68;
t58 = -t22 * t41 + t44 * t34 + r_base(1);
t57 = pkin(2) * t70 + t45 * qJ(3) + t60;
t56 = t24 * t45;
t52 = -pkin(9) - pkin(8);
t39 = qJ(4) + qJ(5);
t37 = cos(t39);
t36 = sin(t39);
t33 = pkin(4) * t50 + pkin(3);
t21 = t59 * t45;
t20 = t59 * t42;
t19 = t40 * t70 - t42 * t68;
t12 = -t21 * t41 + t24 * t44;
t11 = -t41 * t56 - t44 * t59;
t10 = t21 * t44 + t24 * t41;
t9 = -t41 * t59 + t44 * t56;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t62 - mrSges(2,3) - m(3) * t60 - (t48 * mrSges(3,1) + t85) * t42 - t67 * mrSges(5,1) + t84 * t57 + t83 * (pkin(4) * t67 - t19 * t52 + t20 * t33 + t57) + t77 * t20 + (-t82 - t86) * t45 + t79 * (t20 * t36 - t45 * t37) + t76 * (t20 * t37 + t36 * t45) - t74 * t19) * g(3) + (t63 * mrSges(5,1) - mrSges(1,2) + t84 * t61 + t77 * t10 - t75 * t44 + t78 * t41 + t81 * r_base(2) + t80 * t72 + t83 * (-pkin(4) * t63 + t10 * t33 + t9 * t52 + t61) + t79 * (t10 * t36 + t37 * t72) + t76 * (t10 * t37 - t36 * t72) + t74 * t9) * g(2) + (-t64 * mrSges(5,1) - mrSges(1,1) + t84 * t58 + t77 * t12 + t75 * t41 + t78 * t44 + t81 * r_base(1) - t80 * t73 + t83 * (pkin(4) * t64 + t11 * t52 + t12 * t33 + t58) + t79 * (t12 * t36 - t37 * t73) + t76 * (t12 * t37 + t36 * t73) + t74 * t11) * g(1);
U  = t1;
