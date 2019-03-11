% Calculate potential energy for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_energypot_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:03:47
% EndTime: 2019-03-09 01:03:48
% DurationCPUTime: 0.85s
% Computational Cost: add. (584->112), mult. (1402->139), div. (0->0), fcn. (1753->16), ass. (0->56)
t42 = sin(pkin(7));
t45 = cos(pkin(7));
t46 = cos(pkin(6));
t43 = sin(pkin(6));
t52 = cos(qJ(2));
t79 = t43 * t52;
t64 = -t42 * t79 + t45 * t46;
t41 = sin(pkin(13));
t44 = cos(pkin(13));
t50 = sin(qJ(2));
t76 = t46 * t52;
t25 = -t41 * t76 - t44 * t50;
t81 = t43 * t45;
t65 = -t25 * t42 + t41 * t81;
t91 = -m(1) - m(2);
t90 = -m(5) - m(6);
t89 = -m(6) * pkin(11) + m(7) * (-pkin(12) - pkin(11)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t40 = qJ(5) + qJ(6);
t35 = sin(t40);
t36 = cos(t40);
t47 = sin(qJ(5));
t51 = cos(qJ(5));
t88 = -m(7) * (pkin(5) * t47 + pkin(10)) - t47 * mrSges(6,1) - t35 * mrSges(7,1) - t51 * mrSges(6,2) - t36 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t87 = -m(6) * pkin(4) - m(7) * (pkin(5) * t51 + pkin(4)) - t51 * mrSges(6,1) - t36 * mrSges(7,1) + t47 * mrSges(6,2) + t35 * mrSges(7,2) - mrSges(5,1);
t86 = cos(qJ(3));
t85 = cos(qJ(4));
t83 = t41 * t43;
t82 = t43 * t44;
t80 = t43 * t50;
t77 = t46 * t50;
t73 = qJ(1) + r_base(3);
t71 = t42 * t86;
t70 = t45 * t86;
t69 = t44 * pkin(1) + pkin(8) * t83 + r_base(1);
t68 = t46 * pkin(8) + t73;
t67 = t43 * t71;
t23 = -t41 * t50 + t44 * t76;
t66 = -t23 * t42 - t44 * t81;
t63 = t41 * pkin(1) - pkin(8) * t82 + r_base(2);
t26 = -t41 * t77 + t44 * t52;
t62 = t26 * pkin(2) + t65 * pkin(9) + t69;
t49 = sin(qJ(3));
t10 = t26 * t86 + (t25 * t45 + t42 * t83) * t49;
t61 = t10 * pkin(3) + t62;
t60 = pkin(2) * t80 + t64 * pkin(9) + t68;
t17 = t46 * t42 * t49 + (t45 * t49 * t52 + t50 * t86) * t43;
t59 = t17 * pkin(3) + t60;
t24 = t41 * t52 + t44 * t77;
t56 = t24 * pkin(2) + pkin(9) * t66 + t63;
t8 = t24 * t86 + (t23 * t45 - t42 * t82) * t49;
t55 = t8 * pkin(3) + t56;
t48 = sin(qJ(4));
t16 = -t46 * t71 + t49 * t80 - t70 * t79;
t9 = -t25 * t70 + t26 * t49 - t41 * t67;
t7 = -t23 * t70 + t24 * t49 + t44 * t67;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t73 - mrSges(2,3) - m(3) * t68 - t46 * mrSges(3,3) - (t50 * mrSges(3,1) + t52 * mrSges(3,2)) * t43 - m(4) * t60 - t17 * mrSges(4,1) - t64 * mrSges(4,3) - m(7) * t59 + t90 * (t16 * pkin(10) + t59) + t87 * (t17 * t85 + t48 * t64) + t88 * t16 + t89 * (t17 * t48 - t64 * t85)) * g(3) + (-m(3) * t63 - m(4) * t56 - m(7) * t55 - t41 * mrSges(2,1) - t24 * mrSges(3,1) - t8 * mrSges(4,1) - t44 * mrSges(2,2) - t23 * mrSges(3,2) + mrSges(3,3) * t82 - t66 * mrSges(4,3) - mrSges(1,2) + t91 * r_base(2) + t90 * (t7 * pkin(10) + t55) + t87 * (t48 * t66 + t8 * t85) + t88 * t7 + t89 * (t48 * t8 - t66 * t85)) * g(2) + (-m(3) * t69 - m(4) * t62 - m(7) * t61 - t44 * mrSges(2,1) - t26 * mrSges(3,1) - t10 * mrSges(4,1) + t41 * mrSges(2,2) - t25 * mrSges(3,2) - mrSges(3,3) * t83 - t65 * mrSges(4,3) - mrSges(1,1) + t91 * r_base(1) + t90 * (t9 * pkin(10) + t61) + t87 * (t10 * t85 + t48 * t65) + t88 * t9 + t89 * (t10 * t48 - t65 * t85)) * g(1);
U  = t1;
