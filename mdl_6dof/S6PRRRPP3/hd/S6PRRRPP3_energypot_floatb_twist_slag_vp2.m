% Calculate potential energy for
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPP3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:21
% EndTime: 2019-03-08 22:54:21
% DurationCPUTime: 0.60s
% Computational Cost: add. (356->102), mult. (753->114), div. (0->0), fcn. (903->10), ass. (0->56)
t75 = -m(1) - m(2);
t74 = -mrSges(4,3) + mrSges(3,2);
t73 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t72 = mrSges(4,2) - mrSges(5,3) - mrSges(6,1) - m(7) * (pkin(5) + pkin(9)) - mrSges(7,1);
t71 = -m(7) * qJ(6) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t36 = sin(pkin(10));
t38 = cos(pkin(10));
t43 = cos(qJ(2));
t41 = sin(qJ(2));
t61 = cos(pkin(6));
t58 = t41 * t61;
t22 = t36 * t43 + t38 * t58;
t40 = sin(qJ(3));
t37 = sin(pkin(6));
t67 = cos(qJ(3));
t59 = t37 * t67;
t10 = t22 * t40 + t38 * t59;
t70 = pkin(9) * t10;
t24 = -t36 * t58 + t38 * t43;
t12 = t24 * t40 - t36 * t59;
t69 = pkin(9) * t12;
t63 = t37 * t41;
t25 = t40 * t63 - t61 * t67;
t68 = t25 * pkin(9);
t66 = t36 * t37;
t65 = t37 * t38;
t64 = t37 * t40;
t62 = t37 * t43;
t60 = qJ(1) + r_base(3);
t57 = t43 * t61;
t56 = t38 * pkin(1) + pkin(7) * t66 + r_base(1);
t55 = t61 * pkin(7) + t60;
t53 = t36 * pkin(1) - pkin(7) * t65 + r_base(2);
t23 = t36 * t57 + t38 * t41;
t52 = t24 * pkin(2) + pkin(8) * t23 + t56;
t13 = t24 * t67 + t36 * t64;
t51 = t13 * pkin(3) + t52;
t50 = pkin(2) * t63 - pkin(8) * t62 + t55;
t26 = t40 * t61 + t41 * t59;
t49 = t26 * pkin(3) + t50;
t21 = t36 * t41 - t38 * t57;
t48 = t22 * pkin(2) + pkin(8) * t21 + t53;
t11 = t22 * t67 - t38 * t64;
t47 = t11 * pkin(3) + t48;
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t5 = t13 * t39 - t23 * t42;
t6 = t13 * t42 + t23 * t39;
t46 = t6 * pkin(4) + qJ(5) * t5 + t51;
t14 = t26 * t39 + t42 * t62;
t15 = t26 * t42 - t39 * t62;
t45 = t15 * pkin(4) + t14 * qJ(5) + t49;
t3 = t11 * t39 - t21 * t42;
t4 = t11 * t42 + t21 * t39;
t44 = t4 * pkin(4) + qJ(5) * t3 + t47;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t60 - mrSges(2,3) - m(3) * t55 - t61 * mrSges(3,3) - (t41 * mrSges(3,1) + t43 * mrSges(3,2)) * t37 - m(4) * t50 - t26 * mrSges(4,1) + mrSges(4,3) * t62 - m(5) * (t49 + t68) - m(6) * (t45 + t68) - m(7) * t45 + t71 * t15 + t73 * t14 + t72 * t25) * g(3) + (-m(5) * (t47 + t70) - m(6) * (t44 + t70) + mrSges(3,3) * t65 - m(3) * t53 - m(4) * t48 - m(7) * t44 - mrSges(1,2) - t11 * mrSges(4,1) - t22 * mrSges(3,1) - t36 * mrSges(2,1) - t38 * mrSges(2,2) + t75 * r_base(2) + t74 * t21 + t71 * t4 + t73 * t3 + t72 * t10) * g(2) + (-m(5) * (t51 + t69) - m(6) * (t46 + t69) - m(3) * t56 - m(4) * t52 - m(7) * t46 - mrSges(3,3) * t66 - mrSges(1,1) - t13 * mrSges(4,1) - t24 * mrSges(3,1) + t36 * mrSges(2,2) - t38 * mrSges(2,1) + t75 * r_base(1) + t74 * t23 + t71 * t6 + t73 * t5 + t72 * t12) * g(1);
U  = t1;
