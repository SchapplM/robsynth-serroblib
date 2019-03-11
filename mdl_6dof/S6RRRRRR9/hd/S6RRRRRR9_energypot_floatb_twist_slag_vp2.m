% Calculate potential energy for
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRRR9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRRR9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_energypot_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:36
% EndTime: 2019-03-10 05:16:37
% DurationCPUTime: 0.85s
% Computational Cost: add. (584->112), mult. (1402->136), div. (0->0), fcn. (1753->16), ass. (0->57)
t44 = cos(pkin(6));
t49 = sin(qJ(1));
t51 = cos(qJ(2));
t77 = t49 * t51;
t48 = sin(qJ(2));
t52 = cos(qJ(1));
t78 = t48 * t52;
t25 = -t44 * t77 - t78;
t41 = sin(pkin(7));
t43 = cos(pkin(7));
t42 = sin(pkin(6));
t83 = t42 * t49;
t65 = -t25 * t41 + t43 * t83;
t82 = t42 * t51;
t64 = -t41 * t82 + t43 * t44;
t92 = -m(1) - m(2);
t91 = -m(5) - m(6);
t90 = -m(6) * pkin(12) + m(7) * (-pkin(13) - pkin(12)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t40 = qJ(5) + qJ(6);
t35 = sin(t40);
t36 = cos(t40);
t45 = sin(qJ(5));
t50 = cos(qJ(5));
t89 = -m(7) * (pkin(5) * t45 + pkin(11)) - t45 * mrSges(6,1) - t35 * mrSges(7,1) - t50 * mrSges(6,2) - t36 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t88 = -m(6) * pkin(4) - m(7) * (pkin(5) * t50 + pkin(4)) - t50 * mrSges(6,1) - t36 * mrSges(7,1) + t45 * mrSges(6,2) + t35 * mrSges(7,2) - mrSges(5,1);
t87 = cos(qJ(3));
t86 = cos(qJ(4));
t84 = t42 * t48;
t81 = t42 * t52;
t79 = t48 * t49;
t76 = t51 * t52;
t75 = pkin(8) + r_base(3);
t71 = t41 * t87;
t70 = t43 * t87;
t69 = t44 * pkin(9) + t75;
t68 = t52 * pkin(1) + pkin(9) * t83 + r_base(1);
t67 = t42 * t71;
t23 = t44 * t76 - t79;
t66 = -t23 * t41 - t43 * t81;
t63 = t49 * pkin(1) - pkin(9) * t81 + r_base(2);
t26 = -t44 * t79 + t76;
t62 = t26 * pkin(2) + t65 * pkin(10) + t68;
t61 = pkin(2) * t84 + t64 * pkin(10) + t69;
t47 = sin(qJ(3));
t12 = t26 * t87 + (t25 * t43 + t41 * t83) * t47;
t60 = t12 * pkin(3) + t62;
t17 = t44 * t41 * t47 + (t43 * t47 * t51 + t48 * t87) * t42;
t59 = t17 * pkin(3) + t61;
t24 = t44 * t78 + t77;
t56 = t24 * pkin(2) + pkin(10) * t66 + t63;
t10 = t24 * t87 + (t23 * t43 - t41 * t81) * t47;
t55 = t10 * pkin(3) + t56;
t46 = sin(qJ(4));
t16 = -t44 * t71 + t47 * t84 - t70 * t82;
t11 = -t25 * t70 + t26 * t47 - t49 * t67;
t9 = -t23 * t70 + t24 * t47 + t52 * t67;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t75 - mrSges(2,3) - m(3) * t69 - t44 * mrSges(3,3) - (t48 * mrSges(3,1) + t51 * mrSges(3,2)) * t42 - m(4) * t61 - t17 * mrSges(4,1) - t64 * mrSges(4,3) - m(7) * t59 + t91 * (t16 * pkin(11) + t59) + t88 * (t17 * t86 + t46 * t64) + t89 * t16 + t90 * (t17 * t46 - t64 * t86)) * g(3) + (-m(3) * t63 - m(4) * t56 - m(7) * t55 - t49 * mrSges(2,1) - t24 * mrSges(3,1) - t10 * mrSges(4,1) - t52 * mrSges(2,2) - t23 * mrSges(3,2) + mrSges(3,3) * t81 - t66 * mrSges(4,3) - mrSges(1,2) + t92 * r_base(2) + t91 * (t9 * pkin(11) + t55) + t88 * (t10 * t86 + t46 * t66) + t89 * t9 + t90 * (t10 * t46 - t66 * t86)) * g(2) + (-m(3) * t68 - m(4) * t62 - m(7) * t60 - t52 * mrSges(2,1) - t26 * mrSges(3,1) - t12 * mrSges(4,1) + t49 * mrSges(2,2) - t25 * mrSges(3,2) - mrSges(3,3) * t83 - t65 * mrSges(4,3) - mrSges(1,1) + t92 * r_base(1) + t91 * (t11 * pkin(11) + t60) + t88 * (t12 * t86 + t46 * t65) + t89 * t11 + t90 * (t12 * t46 - t65 * t86)) * g(1);
U  = t1;
