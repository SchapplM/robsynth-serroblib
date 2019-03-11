% Calculate potential energy for
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:12:56
% EndTime: 2019-03-09 09:12:57
% DurationCPUTime: 0.67s
% Computational Cost: add. (229->86), mult. (304->77), div. (0->0), fcn. (293->10), ass. (0->42)
t72 = mrSges(3,2) - mrSges(4,3);
t71 = -m(5) * pkin(3) - mrSges(3,1) - mrSges(4,1);
t26 = pkin(10) + qJ(5);
t20 = sin(t26);
t62 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t70 = t62 * t20;
t32 = sin(qJ(2));
t35 = cos(qJ(2));
t29 = cos(pkin(10));
t28 = sin(pkin(10));
t56 = t28 * t32;
t43 = t29 * t35 + t56;
t55 = t28 * t35;
t44 = t29 * t32 - t55;
t21 = cos(t26);
t5 = t20 * t32 + t21 * t35;
t57 = t21 * t32;
t31 = sin(qJ(6));
t34 = cos(qJ(6));
t61 = -m(7) * pkin(5) - t34 * mrSges(7,1) + t31 * mrSges(7,2) - mrSges(6,1);
t65 = -m(6) - m(7);
t69 = t65 * pkin(4) * t56 - t43 * mrSges(5,1) - t44 * mrSges(5,2) + t72 * t32 + t71 * t35 + t5 * t61 - t62 * t57 - mrSges(2,1);
t67 = -m(1) - m(2);
t66 = -m(4) - m(5);
t33 = sin(qJ(1));
t52 = qJ(3) * t32;
t54 = t33 * t35;
t64 = pkin(2) * t54 + t33 * t52;
t59 = -t31 * mrSges(7,1) - t34 * mrSges(7,2) - mrSges(2,2) + mrSges(4,2) + mrSges(3,3) - mrSges(5,3) - mrSges(6,3);
t36 = cos(qJ(1));
t53 = t35 * t36;
t27 = pkin(6) + r_base(3);
t50 = t33 * pkin(1) + r_base(2);
t49 = t32 * pkin(2) + t27;
t48 = t36 * pkin(1) + t33 * pkin(7) + r_base(1);
t47 = t50 + t64;
t42 = -pkin(7) * t36 + t50;
t41 = pkin(2) * t53 + t36 * t52 + t48;
t30 = -pkin(8) - qJ(4);
t18 = pkin(4) * t29 + pkin(3);
t13 = t32 * t18;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - t44 * mrSges(5,1) + t43 * mrSges(5,2) - m(6) * (t13 + t49) - m(7) * (-pkin(4) * t55 + t13) + t61 * (-t20 * t35 + t57) + t62 * t5 + (-m(7) + t66) * (-qJ(3) * t35 + t49) + (-m(6) * (-pkin(4) * t28 - qJ(3)) - t72) * t35 + t71 * t32 + (-m(2) - m(3)) * t27) * g(3) + (-mrSges(1,2) - m(3) * t42 - m(4) * (t42 + t64) - m(5) * t47 + t67 * r_base(2) + t65 * (t18 * t54 + t47) + t54 * t70 + (-m(5) * (-pkin(7) + qJ(4)) + t65 * (-pkin(7) - t30) + t59) * t36 + t69 * t33) * g(2) + (-m(3) * t48 - mrSges(1,1) + t67 * r_base(1) + t66 * t41 + t65 * (t18 * t53 + t33 * t30 + t41) + t53 * t70 + (m(5) * qJ(4) - t59) * t33 + t69 * t36) * g(1);
U  = t1;
