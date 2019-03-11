% Calculate potential energy for
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRP4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:01
% EndTime: 2019-03-09 00:13:01
% DurationCPUTime: 0.71s
% Computational Cost: add. (367->107), mult. (674->124), div. (0->0), fcn. (793->12), ass. (0->52)
t80 = -m(1) - m(2);
t79 = -m(4) - m(5);
t78 = -m(6) - m(7);
t77 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t76 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t48 = cos(qJ(4));
t75 = -t48 * mrSges(5,2) + mrSges(3,2) - mrSges(4,3);
t74 = -m(5) * pkin(9) + mrSges(4,2) - mrSges(7,2) - mrSges(5,3) - mrSges(6,3);
t45 = sin(qJ(4));
t73 = -m(5) * pkin(3) - t48 * mrSges(5,1) + t45 * mrSges(5,2) - mrSges(4,1);
t41 = sin(pkin(11));
t43 = cos(pkin(11));
t47 = sin(qJ(2));
t44 = cos(pkin(6));
t50 = cos(qJ(2));
t63 = t44 * t50;
t23 = t41 * t47 - t43 * t63;
t72 = t23 * t45;
t25 = t41 * t63 + t43 * t47;
t71 = t25 * t45;
t42 = sin(pkin(6));
t70 = t41 * t42;
t69 = t42 * t43;
t46 = sin(qJ(3));
t68 = t42 * t46;
t67 = t42 * t47;
t49 = cos(qJ(3));
t66 = t42 * t49;
t65 = t42 * t50;
t64 = t44 * t47;
t62 = qJ(1) + r_base(3);
t61 = t43 * pkin(1) + pkin(7) * t70 + r_base(1);
t60 = t44 * pkin(7) + t62;
t59 = pkin(2) * t67 + t60;
t58 = t41 * pkin(1) - pkin(7) * t69 + r_base(2);
t26 = -t41 * t64 + t43 * t50;
t57 = t26 * pkin(2) + pkin(8) * t25 + t61;
t56 = -pkin(8) * t65 + t59;
t24 = t41 * t50 + t43 * t64;
t55 = t24 * pkin(2) + pkin(8) * t23 + t58;
t51 = -pkin(10) - pkin(9);
t40 = qJ(4) + qJ(5);
t36 = cos(t40);
t35 = sin(t40);
t34 = pkin(4) * t48 + pkin(3);
t28 = t44 * t46 + t47 * t66;
t27 = -t44 * t49 + t46 * t67;
t14 = t26 * t49 + t41 * t68;
t13 = t26 * t46 - t41 * t66;
t12 = t24 * t49 - t43 * t68;
t11 = t24 * t46 + t43 * t66;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t62 - mrSges(2,3) - m(3) * t60 - t44 * mrSges(3,3) - (t47 * mrSges(3,1) + t50 * mrSges(3,2)) * t42 - m(4) * t56 - t28 * mrSges(4,1) + mrSges(4,3) * t65 - m(5) * (pkin(3) * t28 + t56) - (t28 * t48 - t45 * t65) * mrSges(5,1) - (-t28 * t45 - t48 * t65) * mrSges(5,2) + t78 * (t28 * t34 - t27 * t51 + (-pkin(4) * t45 - pkin(8)) * t65 + t59) + t77 * (t28 * t36 - t35 * t65) + t76 * (t28 * t35 + t36 * t65) + t74 * t27) * g(3) + (-m(3) * t58 - t41 * mrSges(2,1) - t24 * mrSges(3,1) - t72 * mrSges(5,1) - t43 * mrSges(2,2) + mrSges(3,3) * t69 - mrSges(1,2) + t80 * r_base(2) + t79 * t55 + t78 * (pkin(4) * t72 - t11 * t51 + t12 * t34 + t55) + t73 * t12 + t75 * t23 + t77 * (t12 * t36 + t23 * t35) + t76 * (t12 * t35 - t23 * t36) + t74 * t11) * g(2) + (-m(3) * t61 - t43 * mrSges(2,1) - t26 * mrSges(3,1) - t71 * mrSges(5,1) + t41 * mrSges(2,2) - mrSges(3,3) * t70 - mrSges(1,1) + t80 * r_base(1) + t79 * t57 + t78 * (pkin(4) * t71 - t13 * t51 + t14 * t34 + t57) + t77 * (t14 * t36 + t25 * t35) + t76 * (t14 * t35 - t25 * t36) + t73 * t14 + t75 * t25 + t74 * t13) * g(1);
U  = t1;
