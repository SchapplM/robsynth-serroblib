% Calculate potential energy for
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRP10_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:22
% EndTime: 2019-03-09 17:29:22
% DurationCPUTime: 0.71s
% Computational Cost: add. (367->107), mult. (674->121), div. (0->0), fcn. (793->12), ass. (0->53)
t81 = -m(1) - m(2);
t80 = -m(4) - m(5);
t79 = -m(6) - m(7);
t78 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t77 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t43 = cos(pkin(11));
t76 = -t43 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t75 = -m(5) * qJ(4) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t41 = sin(pkin(11));
t74 = -m(5) * pkin(3) - t43 * mrSges(5,1) + t41 * mrSges(5,2) - mrSges(4,1);
t44 = cos(pkin(6));
t50 = cos(qJ(2));
t51 = cos(qJ(1));
t63 = t50 * t51;
t47 = sin(qJ(2));
t48 = sin(qJ(1));
t66 = t47 * t48;
t25 = -t44 * t63 + t66;
t73 = t25 * t41;
t64 = t48 * t50;
t65 = t47 * t51;
t27 = t44 * t64 + t65;
t72 = t27 * t41;
t42 = sin(pkin(6));
t71 = t42 * t47;
t70 = t42 * t48;
t49 = cos(qJ(3));
t69 = t42 * t49;
t68 = t42 * t50;
t67 = t42 * t51;
t62 = pkin(7) + r_base(3);
t61 = t44 * pkin(8) + t62;
t60 = t51 * pkin(1) + pkin(8) * t70 + r_base(1);
t59 = pkin(2) * t71 + t61;
t58 = t48 * pkin(1) - pkin(8) * t67 + r_base(2);
t28 = -t44 * t66 + t63;
t57 = t28 * pkin(2) + pkin(9) * t27 + t60;
t56 = -pkin(9) * t68 + t59;
t26 = t44 * t65 + t64;
t55 = t26 * pkin(2) + t25 * pkin(9) + t58;
t46 = sin(qJ(3));
t45 = -pkin(10) - qJ(4);
t40 = pkin(11) + qJ(5);
t36 = cos(t40);
t35 = sin(t40);
t34 = pkin(4) * t43 + pkin(3);
t24 = t44 * t46 + t47 * t69;
t23 = -t44 * t49 + t46 * t71;
t14 = t28 * t49 + t46 * t70;
t13 = t28 * t46 - t48 * t69;
t12 = t26 * t49 - t46 * t67;
t11 = t26 * t46 + t49 * t67;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t62 - mrSges(2,3) - m(3) * t61 - t44 * mrSges(3,3) - (t47 * mrSges(3,1) + t50 * mrSges(3,2)) * t42 - m(4) * t56 - t24 * mrSges(4,1) + mrSges(4,3) * t68 - m(5) * (pkin(3) * t24 + t56) - (t24 * t43 - t41 * t68) * mrSges(5,1) - (-t24 * t41 - t43 * t68) * mrSges(5,2) + t79 * (t24 * t34 - t23 * t45 + (-pkin(4) * t41 - pkin(9)) * t68 + t59) + t78 * (t24 * t36 - t35 * t68) + t77 * (t24 * t35 + t36 * t68) + t75 * t23) * g(3) + (-m(3) * t58 - t48 * mrSges(2,1) - t26 * mrSges(3,1) - t73 * mrSges(5,1) - t51 * mrSges(2,2) + mrSges(3,3) * t67 - mrSges(1,2) + t81 * r_base(2) + t80 * t55 + t79 * (pkin(4) * t73 - t11 * t45 + t12 * t34 + t55) + t74 * t12 + t76 * t25 + t78 * (t12 * t36 + t25 * t35) + t77 * (t12 * t35 - t25 * t36) + t75 * t11) * g(2) + (-m(3) * t60 - t51 * mrSges(2,1) - t28 * mrSges(3,1) - t72 * mrSges(5,1) + t48 * mrSges(2,2) - mrSges(3,3) * t70 - mrSges(1,1) + t81 * r_base(1) + t80 * t57 + t79 * (pkin(4) * t72 - t13 * t45 + t14 * t34 + t57) + t78 * (t14 * t36 + t27 * t35) + t77 * (t14 * t35 - t27 * t36) + t74 * t14 + t76 * t27 + t75 * t13) * g(1);
U  = t1;
