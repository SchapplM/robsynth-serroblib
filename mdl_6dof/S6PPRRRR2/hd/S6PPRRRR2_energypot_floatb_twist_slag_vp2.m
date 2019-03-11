% Calculate potential energy for
% S6PPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PPRRRR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PPRRRR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_energypot_floatb_twist_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:09
% EndTime: 2019-03-08 19:03:10
% DurationCPUTime: 0.86s
% Computational Cost: add. (584->112), mult. (1402->140), div. (0->0), fcn. (1753->16), ass. (0->57)
t43 = sin(pkin(7));
t47 = cos(pkin(7));
t48 = cos(pkin(6));
t44 = sin(pkin(6));
t45 = cos(pkin(13));
t81 = t44 * t45;
t64 = -t43 * t81 + t47 * t48;
t41 = sin(pkin(13));
t46 = cos(pkin(12));
t42 = sin(pkin(12));
t82 = t42 * t48;
t25 = -t41 * t46 - t45 * t82;
t79 = t44 * t47;
t65 = -t25 * t43 + t42 * t79;
t92 = -m(1) - m(2);
t91 = -m(5) - m(6);
t90 = -m(6) * pkin(10) + m(7) * (-pkin(11) - pkin(10)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t40 = qJ(5) + qJ(6);
t36 = sin(t40);
t37 = cos(t40);
t49 = sin(qJ(5));
t52 = cos(qJ(5));
t89 = -m(7) * (pkin(5) * t49 + pkin(9)) - t49 * mrSges(6,1) - t36 * mrSges(7,1) - t52 * mrSges(6,2) - t37 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t88 = -m(6) * pkin(4) - m(7) * (pkin(5) * t52 + pkin(4)) - t52 * mrSges(6,1) - t37 * mrSges(7,1) + t49 * mrSges(6,2) + t36 * mrSges(7,2) - mrSges(5,1);
t87 = cos(qJ(3));
t86 = cos(qJ(4));
t84 = t41 * t44;
t83 = t42 * t44;
t80 = t44 * t46;
t78 = t46 * t48;
t76 = qJ(2) * t44;
t73 = qJ(1) + r_base(3);
t71 = t43 * t87;
t70 = t47 * t87;
t69 = t46 * pkin(1) + t42 * t76 + r_base(1);
t68 = t48 * qJ(2) + t73;
t67 = t44 * t71;
t23 = -t41 * t42 + t45 * t78;
t66 = -t23 * t43 - t46 * t79;
t63 = t42 * pkin(1) - t46 * t76 + r_base(2);
t26 = -t41 * t82 + t45 * t46;
t62 = t26 * pkin(2) + t65 * pkin(8) + t69;
t51 = sin(qJ(3));
t10 = t26 * t87 + (t25 * t47 + t43 * t83) * t51;
t61 = t10 * pkin(3) + t62;
t60 = pkin(2) * t84 + t64 * pkin(8) + t68;
t17 = t48 * t43 * t51 + (t45 * t47 * t51 + t41 * t87) * t44;
t59 = t17 * pkin(3) + t60;
t24 = t41 * t78 + t42 * t45;
t56 = t24 * pkin(2) + pkin(8) * t66 + t63;
t8 = t24 * t87 + (t23 * t47 - t43 * t80) * t51;
t55 = t8 * pkin(3) + t56;
t50 = sin(qJ(4));
t16 = -t48 * t71 + t51 * t84 - t70 * t81;
t9 = -t25 * t70 + t26 * t51 - t42 * t67;
t7 = -t23 * t70 + t24 * t51 + t46 * t67;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t73 - mrSges(2,3) - m(3) * t68 - t48 * mrSges(3,3) - (t41 * mrSges(3,1) + t45 * mrSges(3,2)) * t44 - m(4) * t60 - t17 * mrSges(4,1) - t64 * mrSges(4,3) - m(7) * t59 + t91 * (t16 * pkin(9) + t59) + t88 * (t17 * t86 + t50 * t64) + t89 * t16 + t90 * (t17 * t50 - t64 * t86)) * g(3) + (-m(3) * t63 - m(4) * t56 - m(7) * t55 - t42 * mrSges(2,1) - t24 * mrSges(3,1) - t8 * mrSges(4,1) - t46 * mrSges(2,2) - t23 * mrSges(3,2) + mrSges(3,3) * t80 - t66 * mrSges(4,3) - mrSges(1,2) + t92 * r_base(2) + t91 * (t7 * pkin(9) + t55) + t88 * (t50 * t66 + t8 * t86) + t89 * t7 + t90 * (t50 * t8 - t66 * t86)) * g(2) + (-m(3) * t69 - m(4) * t62 - m(7) * t61 - t46 * mrSges(2,1) - t26 * mrSges(3,1) - t10 * mrSges(4,1) + t42 * mrSges(2,2) - t25 * mrSges(3,2) - mrSges(3,3) * t83 - t65 * mrSges(4,3) - mrSges(1,1) + t92 * r_base(1) + t91 * (t9 * pkin(9) + t61) + t88 * (t10 * t86 + t50 * t65) + t89 * t9 + t90 * (t10 * t50 - t65 * t86)) * g(1);
U  = t1;
