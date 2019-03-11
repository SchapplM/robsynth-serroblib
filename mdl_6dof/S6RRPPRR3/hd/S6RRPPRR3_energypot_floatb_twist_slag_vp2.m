% Calculate potential energy for
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPPRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:24
% EndTime: 2019-03-09 08:54:25
% DurationCPUTime: 0.82s
% Computational Cost: add. (404->102), mult. (738->120), div. (0->0), fcn. (880->14), ass. (0->56)
t83 = -m(4) - m(5);
t82 = -m(6) - m(7);
t81 = -m(2) - m(3) - m(1);
t80 = -m(3) * pkin(1) - mrSges(2,1);
t79 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t43 = cos(pkin(12));
t78 = -t43 * mrSges(5,2) - mrSges(3,3) - mrSges(4,3);
t77 = m(3) * pkin(8) - t78;
t40 = sin(pkin(12));
t76 = -m(5) * pkin(3) - t43 * mrSges(5,1) + t40 * mrSges(5,2) - mrSges(4,1);
t47 = sin(qJ(6));
t50 = cos(qJ(6));
t75 = -m(7) * pkin(5) - t50 * mrSges(7,1) + t47 * mrSges(7,2) - mrSges(6,1);
t74 = m(5) * qJ(4) + t47 * mrSges(7,1) + t50 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t45 = cos(pkin(6));
t73 = t40 * t45;
t42 = sin(pkin(6));
t48 = sin(qJ(2));
t72 = t42 * t48;
t49 = sin(qJ(1));
t71 = t42 * t49;
t52 = cos(qJ(1));
t70 = t42 * t52;
t44 = cos(pkin(11));
t51 = cos(qJ(2));
t69 = t44 * t51;
t68 = t48 * t52;
t67 = t49 * t48;
t66 = t49 * t51;
t65 = t51 * t52;
t64 = pkin(7) + r_base(3);
t63 = t40 * t71;
t62 = t40 * t70;
t61 = t45 * pkin(8) + t64;
t22 = pkin(2) * t45 * t48 + (-pkin(8) - qJ(3)) * t42;
t34 = pkin(2) * t51 + pkin(1);
t60 = t52 * t22 + t49 * t34 + r_base(2);
t41 = sin(pkin(11));
t59 = t41 * t51 + t44 * t48;
t24 = -t41 * t48 + t69;
t58 = -t22 * t49 + t52 * t34 + r_base(1);
t57 = pkin(2) * t72 + t45 * qJ(3) + t61;
t56 = t24 * t45;
t46 = -pkin(9) - qJ(4);
t39 = pkin(12) + qJ(5);
t36 = cos(t39);
t35 = sin(t39);
t33 = pkin(4) * t43 + pkin(3);
t21 = t59 * t45;
t20 = t59 * t42;
t19 = t41 * t72 - t42 * t69;
t12 = -t49 * t21 + t24 * t52;
t11 = -t49 * t56 - t52 * t59;
t10 = t21 * t52 + t49 * t24;
t9 = -t49 * t59 + t52 * t56;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t64 - mrSges(2,3) - m(3) * t61 - (t48 * mrSges(3,1) + t51 * mrSges(3,2)) * t42 - t73 * mrSges(5,1) + t83 * t57 + t82 * (pkin(4) * t73 - t19 * t46 + t20 * t33 + t57) + t76 * t20 + t78 * t45 + t79 * (t20 * t35 - t45 * t36) + t75 * (t20 * t36 + t35 * t45) - t74 * t19) * g(3) + (-mrSges(1,2) + t62 * mrSges(5,1) - (t45 * t65 - t67) * mrSges(3,2) - (t45 * t68 + t66) * mrSges(3,1) - t52 * mrSges(2,2) + t80 * t49 + t81 * r_base(2) + t83 * t60 + t76 * t10 + t77 * t70 + t82 * (-pkin(4) * t62 + t10 * t33 + t9 * t46 + t60) + t79 * (t10 * t35 + t36 * t70) + t75 * (t10 * t36 - t35 * t70) + t74 * t9) * g(2) + (-mrSges(1,1) - t63 * mrSges(5,1) - (-t45 * t67 + t65) * mrSges(3,1) - (-t45 * t66 - t68) * mrSges(3,2) + t49 * mrSges(2,2) + t80 * t52 + t81 * r_base(1) + t83 * t58 + t76 * t12 - t77 * t71 + t82 * (pkin(4) * t63 + t11 * t46 + t12 * t33 + t58) + t79 * (t12 * t35 - t36 * t71) + t75 * (t12 * t36 + t35 * t71) + t74 * t11) * g(1);
U  = t1;
