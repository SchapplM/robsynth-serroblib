% Calculate potential energy for
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:32
% EndTime: 2019-03-08 19:24:33
% DurationCPUTime: 0.85s
% Computational Cost: add. (404->98), mult. (738->113), div. (0->0), fcn. (880->14), ass. (0->55)
t51 = cos(qJ(4));
t87 = t51 * mrSges(5,2);
t52 = cos(qJ(2));
t86 = t52 * mrSges(3,2);
t85 = -m(4) - m(5);
t84 = -m(6) - m(7);
t83 = mrSges(3,3) + mrSges(4,3);
t40 = sin(pkin(11));
t43 = cos(pkin(11));
t49 = sin(qJ(2));
t82 = t49 * t40 - t43 * t52;
t81 = -m(3) - m(1) - m(2);
t80 = m(3) * pkin(7) + t83;
t79 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t78 = -m(3) * pkin(1) - t52 * mrSges(3,1) + t49 * mrSges(3,2) - mrSges(2,1);
t48 = sin(qJ(4));
t77 = -m(5) * pkin(3) - t51 * mrSges(5,1) + t48 * mrSges(5,2) - mrSges(4,1);
t47 = sin(qJ(6));
t50 = cos(qJ(6));
t76 = -m(7) * pkin(5) - t50 * mrSges(7,1) + t47 * mrSges(7,2) - mrSges(6,1);
t42 = sin(pkin(6));
t45 = cos(pkin(6));
t67 = t45 * t49;
t75 = t67 * mrSges(3,1) - t42 * t87 + t45 * t86 + mrSges(2,2);
t74 = m(5) * pkin(8) + t47 * mrSges(7,1) + t50 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t41 = sin(pkin(10));
t73 = t41 * t42;
t44 = cos(pkin(10));
t72 = t42 * t44;
t71 = t42 * t48;
t68 = t45 * t48;
t64 = t41 * t71;
t63 = t44 * t71;
t62 = qJ(1) + r_base(3);
t22 = pkin(2) * t67 + (-pkin(7) - qJ(3)) * t42;
t34 = pkin(2) * t52 + pkin(1);
t61 = t44 * t22 + t41 * t34 + r_base(2);
t60 = t45 * pkin(7) + t62;
t59 = t40 * t52 + t49 * t43;
t58 = -t22 * t41 + t44 * t34 + r_base(1);
t57 = t42 * t49 * pkin(2) + t45 * qJ(3) + t60;
t56 = t82 * t45;
t46 = -qJ(5) - pkin(8);
t39 = qJ(4) + pkin(12);
t36 = cos(t39);
t35 = sin(t39);
t33 = pkin(4) * t51 + pkin(3);
t21 = t59 * t45;
t20 = t59 * t42;
t19 = t82 * t42;
t12 = -t21 * t41 - t44 * t82;
t11 = t41 * t56 - t44 * t59;
t10 = t21 * t44 - t41 * t82;
t9 = -t41 * t59 - t44 * t56;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t62 - mrSges(2,3) - m(3) * t60 - (t49 * mrSges(3,1) + t86) * t42 - t68 * mrSges(5,1) + t85 * t57 + t84 * (pkin(4) * t68 - t19 * t46 + t20 * t33 + t57) + t77 * t20 + (-t83 - t87) * t45 + t79 * (t20 * t35 - t45 * t36) + t76 * (t20 * t36 + t35 * t45) - t74 * t19) * g(3) + (t63 * mrSges(5,1) - mrSges(1,2) + t85 * t61 + t77 * t10 - t75 * t44 + t78 * t41 + t81 * r_base(2) + t80 * t72 + t84 * (-pkin(4) * t63 + t10 * t33 + t9 * t46 + t61) + t79 * (t10 * t35 + t36 * t72) + t76 * (t10 * t36 - t35 * t72) + t74 * t9) * g(2) + (-t64 * mrSges(5,1) - mrSges(1,1) + t85 * t58 + t77 * t12 + t75 * t41 + t78 * t44 + t81 * r_base(1) - t80 * t73 + t84 * (pkin(4) * t64 + t11 * t46 + t12 * t33 + t58) + t79 * (t12 * t35 - t36 * t73) + t76 * (t12 * t36 + t35 * t73) + t74 * t11) * g(1);
U  = t1;
