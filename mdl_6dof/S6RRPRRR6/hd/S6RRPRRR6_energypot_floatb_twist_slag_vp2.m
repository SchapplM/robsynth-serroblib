% Calculate potential energy for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:50:22
% EndTime: 2019-03-09 13:50:23
% DurationCPUTime: 0.67s
% Computational Cost: add. (229->86), mult. (304->77), div. (0->0), fcn. (293->10), ass. (0->42)
t72 = mrSges(3,2) - mrSges(4,3);
t71 = -m(5) * pkin(3) - mrSges(3,1) - mrSges(4,1);
t27 = qJ(4) + qJ(5);
t20 = sin(t27);
t62 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t70 = t62 * t20;
t30 = sin(qJ(2));
t34 = cos(qJ(2));
t33 = cos(qJ(4));
t29 = sin(qJ(4));
t56 = t29 * t30;
t43 = t33 * t34 + t56;
t55 = t29 * t34;
t44 = t30 * t33 - t55;
t21 = cos(t27);
t5 = t20 * t30 + t21 * t34;
t57 = t21 * t30;
t28 = sin(qJ(6));
t32 = cos(qJ(6));
t61 = -m(7) * pkin(5) - mrSges(7,1) * t32 + mrSges(7,2) * t28 - mrSges(6,1);
t65 = -m(6) - m(7);
t69 = pkin(4) * t56 * t65 - t43 * mrSges(5,1) - t44 * mrSges(5,2) + t30 * t72 + t34 * t71 + t5 * t61 - t57 * t62 - mrSges(2,1);
t67 = -m(1) - m(2);
t66 = -m(4) - m(5);
t31 = sin(qJ(1));
t52 = qJ(3) * t30;
t54 = t31 * t34;
t64 = pkin(2) * t54 + t31 * t52;
t59 = -t28 * mrSges(7,1) - t32 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3) - mrSges(6,3);
t35 = cos(qJ(1));
t53 = t34 * t35;
t26 = pkin(6) + r_base(3);
t50 = pkin(1) * t31 + r_base(2);
t49 = pkin(2) * t30 + t26;
t48 = pkin(1) * t35 + pkin(7) * t31 + r_base(1);
t47 = t50 + t64;
t42 = -pkin(7) * t35 + t50;
t41 = pkin(2) * t53 + t35 * t52 + t48;
t36 = -pkin(9) - pkin(8);
t18 = pkin(4) * t33 + pkin(3);
t13 = t30 * t18;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t44 * mrSges(5,1) + t43 * mrSges(5,2) - m(6) * (t13 + t49) - m(7) * (-pkin(4) * t55 + t13) + t61 * (-t20 * t34 + t57) + t62 * t5 + (-m(7) + t66) * (-qJ(3) * t34 + t49) + (-m(6) * (-pkin(4) * t29 - qJ(3)) - t72) * t34 + t71 * t30 + (-m(2) - m(3)) * t26) * g(3) + (-mrSges(1,2) - m(3) * t42 - m(4) * (t42 + t64) - m(5) * t47 + t67 * r_base(2) + t65 * (t18 * t54 + t47) + t54 * t70 + (-m(5) * (-pkin(7) + pkin(8)) + t65 * (-pkin(7) - t36) + t59) * t35 + t69 * t31) * g(2) + (-m(3) * t48 - mrSges(1,1) + t67 * r_base(1) + t66 * t41 + t65 * (t18 * t53 + t31 * t36 + t41) + t53 * t70 + (m(5) * pkin(8) - t59) * t31 + t69 * t35) * g(1);
U  = t1;
