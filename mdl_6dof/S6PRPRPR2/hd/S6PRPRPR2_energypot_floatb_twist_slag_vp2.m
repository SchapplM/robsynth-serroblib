% Calculate potential energy for
% S6PRPRPR2
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:25
% EndTime: 2019-03-08 19:29:26
% DurationCPUTime: 0.83s
% Computational Cost: add. (421->91), mult. (867->97), div. (0->0), fcn. (1059->14), ass. (0->47)
t46 = cos(qJ(2));
t80 = t46 * mrSges(3,2);
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t34 = pkin(12) + qJ(6);
t30 = sin(t34);
t31 = cos(t34);
t35 = sin(pkin(12));
t39 = cos(pkin(12));
t69 = -m(6) * pkin(4) - m(7) * (pkin(5) * t39 + pkin(4)) - t39 * mrSges(6,1) - t31 * mrSges(7,1) + t35 * mrSges(6,2) + t30 * mrSges(7,2) - mrSges(5,1);
t72 = -m(6) * qJ(5) + m(7) * (-pkin(9) - qJ(5)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t79 = t43 * t72 + t45 * t69 - mrSges(4,1);
t38 = sin(pkin(6));
t41 = cos(pkin(6));
t44 = sin(qJ(2));
t64 = t41 * t44;
t76 = -mrSges(3,3) - mrSges(4,3);
t78 = -t64 * mrSges(3,1) - t41 * t80 - mrSges(2,2) + (m(3) * pkin(7) - t43 * t69 + t45 * t72 - t76) * t38;
t77 = -m(5) - m(6);
t36 = sin(pkin(11));
t61 = cos(pkin(11));
t21 = -t44 * t36 + t46 * t61;
t75 = -m(3) - m(1) - m(2);
t71 = -m(3) * pkin(1) - t46 * mrSges(3,1) + t44 * mrSges(3,2) - mrSges(2,1);
t70 = m(7) * (pkin(5) * t35 + pkin(8)) + t35 * mrSges(6,1) + t30 * mrSges(7,1) + t39 * mrSges(6,2) + t31 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3);
t60 = qJ(1) + r_base(3);
t19 = pkin(2) * t64 + (-pkin(7) - qJ(3)) * t38;
t29 = pkin(2) * t46 + pkin(1);
t37 = sin(pkin(10));
t40 = cos(pkin(10));
t57 = t40 * t19 + t37 * t29 + r_base(2);
t56 = t41 * pkin(7) + t60;
t20 = -t46 * t36 - t44 * t61;
t18 = t20 * t41;
t8 = -t18 * t40 + t21 * t37;
t55 = t8 * pkin(3) + t57;
t54 = -t19 * t37 + t40 * t29 + r_base(1);
t53 = t38 * t44 * pkin(2) + t41 * qJ(3) + t56;
t10 = t18 * t37 + t21 * t40;
t52 = t10 * pkin(3) + t54;
t17 = t20 * t38;
t51 = -t17 * pkin(3) + t53;
t48 = t41 * t21;
t16 = t21 * t38;
t9 = t20 * t40 - t37 * t48;
t7 = t37 * t20 + t40 * t48;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t60 - mrSges(2,3) - m(3) * t56 - (t44 * mrSges(3,1) + t80) * t38 - m(4) * t53 + t17 * mrSges(4,1) - m(7) * t51 + t77 * (-pkin(8) * t16 + t51) + t76 * t41 + t69 * (-t17 * t45 + t41 * t43) + t70 * t16 + t72 * (-t17 * t43 - t41 * t45)) * g(3) + (-m(4) * t57 - m(7) * t55 - mrSges(1,2) + t71 * t37 + t75 * r_base(2) + t77 * (-pkin(8) * t7 + t55) + t70 * t7 + t79 * t8 + t78 * t40) * g(2) + (-m(4) * t54 - m(7) * t52 - mrSges(1,1) + t71 * t40 + t75 * r_base(1) + t77 * (-pkin(8) * t9 + t52) + t70 * t9 + t79 * t10 - t78 * t37) * g(1);
U  = t1;
