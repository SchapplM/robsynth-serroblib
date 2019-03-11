% Calculate potential energy for
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:21:08
% EndTime: 2019-03-09 10:21:09
% DurationCPUTime: 0.82s
% Computational Cost: add. (404->102), mult. (738->120), div. (0->0), fcn. (880->14), ass. (0->56)
t83 = -m(4) - m(5);
t82 = -m(6) - m(7);
t81 = -m(3) - m(2) - m(1);
t80 = -m(3) * pkin(1) - mrSges(2,1);
t79 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t50 = cos(qJ(4));
t78 = -t50 * mrSges(5,2) - mrSges(3,3) - mrSges(4,3);
t77 = m(3) * pkin(8) - t78;
t46 = sin(qJ(4));
t76 = -m(5) * pkin(3) - t50 * mrSges(5,1) + t46 * mrSges(5,2) - mrSges(4,1);
t45 = sin(qJ(6));
t49 = cos(qJ(6));
t75 = -m(7) * pkin(5) - t49 * mrSges(7,1) + t45 * mrSges(7,2) - mrSges(6,1);
t74 = m(5) * pkin(9) + t45 * mrSges(7,1) + t49 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t41 = sin(pkin(6));
t47 = sin(qJ(2));
t73 = t41 * t47;
t48 = sin(qJ(1));
t72 = t41 * t48;
t52 = cos(qJ(1));
t71 = t41 * t52;
t42 = cos(pkin(11));
t51 = cos(qJ(2));
t70 = t42 * t51;
t43 = cos(pkin(6));
t69 = t43 * t46;
t68 = t47 * t52;
t67 = t48 * t47;
t66 = t48 * t51;
t65 = t51 * t52;
t64 = pkin(7) + r_base(3);
t63 = t46 * t72;
t62 = t46 * t71;
t61 = t43 * pkin(8) + t64;
t22 = pkin(2) * t43 * t47 + (-pkin(8) - qJ(3)) * t41;
t34 = pkin(2) * t51 + pkin(1);
t60 = t52 * t22 + t48 * t34 + r_base(2);
t40 = sin(pkin(11));
t59 = t40 * t51 + t42 * t47;
t24 = -t40 * t47 + t70;
t58 = -t22 * t48 + t52 * t34 + r_base(1);
t57 = pkin(2) * t73 + t43 * qJ(3) + t61;
t56 = t24 * t43;
t44 = -qJ(5) - pkin(9);
t39 = qJ(4) + pkin(12);
t36 = cos(t39);
t35 = sin(t39);
t33 = pkin(4) * t50 + pkin(3);
t21 = t59 * t43;
t20 = t59 * t41;
t19 = t40 * t73 - t41 * t70;
t12 = -t48 * t21 + t24 * t52;
t11 = -t48 * t56 - t52 * t59;
t10 = t21 * t52 + t48 * t24;
t9 = -t48 * t59 + t52 * t56;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t64 - mrSges(2,3) - m(3) * t61 - (t47 * mrSges(3,1) + t51 * mrSges(3,2)) * t41 - t69 * mrSges(5,1) + t83 * t57 + t82 * (pkin(4) * t69 - t19 * t44 + t20 * t33 + t57) + t76 * t20 + t78 * t43 + t79 * (t20 * t35 - t43 * t36) + t75 * (t20 * t36 + t35 * t43) - t74 * t19) * g(3) + (-(t43 * t65 - t67) * mrSges(3,2) - (t43 * t68 + t66) * mrSges(3,1) + t62 * mrSges(5,1) - mrSges(1,2) - t52 * mrSges(2,2) + t80 * t48 + t81 * r_base(2) + t83 * t60 + t76 * t10 + t77 * t71 + t82 * (-pkin(4) * t62 + t10 * t33 + t9 * t44 + t60) + t79 * (t10 * t35 + t36 * t71) + t75 * (t10 * t36 - t35 * t71) + t74 * t9) * g(2) + (-(-t43 * t66 - t68) * mrSges(3,2) - (-t43 * t67 + t65) * mrSges(3,1) - t63 * mrSges(5,1) - mrSges(1,1) + t48 * mrSges(2,2) + t80 * t52 + t81 * r_base(1) + t83 * t58 + t76 * t12 - t77 * t72 + t82 * (pkin(4) * t63 + t11 * t44 + t12 * t33 + t58) + t79 * (t12 * t35 - t36 * t72) + t75 * (t12 * t36 + t35 * t72) + t74 * t11) * g(1);
U  = t1;
