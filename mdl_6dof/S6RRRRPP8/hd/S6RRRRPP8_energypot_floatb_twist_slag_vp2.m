% Calculate potential energy for
% S6RRRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP8_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:12
% EndTime: 2019-03-09 21:30:12
% DurationCPUTime: 0.60s
% Computational Cost: add. (356->102), mult. (753->113), div. (0->0), fcn. (903->10), ass. (0->55)
t73 = -m(1) - m(2);
t72 = -mrSges(4,3) + mrSges(3,2);
t71 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t70 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2) - m(7) * (pkin(10) - qJ(6)) + mrSges(7,3);
t69 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t38 = sin(qJ(2));
t41 = cos(qJ(2));
t42 = cos(qJ(1));
t39 = sin(qJ(1));
t60 = cos(pkin(6));
t56 = t39 * t60;
t26 = -t38 * t56 + t42 * t41;
t37 = sin(qJ(3));
t35 = sin(pkin(6));
t65 = cos(qJ(3));
t58 = t35 * t65;
t14 = t26 * t37 - t39 * t58;
t68 = pkin(10) * t14;
t64 = t35 * t38;
t21 = t37 * t64 - t60 * t65;
t67 = pkin(10) * t21;
t55 = t42 * t60;
t24 = t38 * t55 + t39 * t41;
t12 = t24 * t37 + t42 * t58;
t66 = t12 * pkin(10);
t63 = t35 * t39;
t62 = t35 * t41;
t61 = t35 * t42;
t59 = pkin(7) + r_base(3);
t57 = t60 * pkin(8) + t59;
t54 = t42 * pkin(1) + pkin(8) * t63 + r_base(1);
t52 = t39 * pkin(1) - pkin(8) * t61 + r_base(2);
t25 = t42 * t38 + t41 * t56;
t51 = t26 * pkin(2) + pkin(9) * t25 + t54;
t50 = pkin(2) * t64 - pkin(9) * t62 + t57;
t15 = t26 * t65 + t37 * t63;
t49 = t15 * pkin(3) + t51;
t22 = t37 * t60 + t38 * t58;
t48 = t22 * pkin(3) + t50;
t23 = t38 * t39 - t41 * t55;
t47 = t24 * pkin(2) + t23 * pkin(9) + t52;
t13 = t24 * t65 - t37 * t61;
t46 = t13 * pkin(3) + t47;
t36 = sin(qJ(4));
t40 = cos(qJ(4));
t5 = t15 * t36 - t25 * t40;
t6 = t15 * t40 + t25 * t36;
t45 = t6 * pkin(4) + qJ(5) * t5 + t49;
t10 = t22 * t36 + t40 * t62;
t11 = t22 * t40 - t36 * t62;
t44 = t11 * pkin(4) + qJ(5) * t10 + t48;
t3 = t13 * t36 - t23 * t40;
t4 = t13 * t40 + t23 * t36;
t43 = t4 * pkin(4) + t3 * qJ(5) + t46;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t59 - mrSges(2,3) - m(3) * t57 - t60 * mrSges(3,3) - (t38 * mrSges(3,1) + t41 * mrSges(3,2)) * t35 - m(4) * t50 - t22 * mrSges(4,1) + mrSges(4,3) * t62 - m(5) * (t48 + t67) - m(6) * (t44 + t67) - m(7) * t44 + t69 * t11 + t71 * t10 + t70 * t21) * g(3) + (-m(5) * (t46 + t66) - m(6) * (t43 + t66) + mrSges(3,3) * t61 - m(3) * t52 - m(4) * t47 - m(7) * t43 - mrSges(1,2) - t13 * mrSges(4,1) - t24 * mrSges(3,1) - t39 * mrSges(2,1) - t42 * mrSges(2,2) + t73 * r_base(2) + t72 * t23 + t69 * t4 + t71 * t3 + t70 * t12) * g(2) + (-m(5) * (t49 + t68) - m(6) * (t45 + t68) - m(3) * t54 - m(7) * t45 - m(4) * t51 - mrSges(1,1) - mrSges(3,3) * t63 - t15 * mrSges(4,1) - t26 * mrSges(3,1) + t39 * mrSges(2,2) - t42 * mrSges(2,1) + t73 * r_base(1) + t72 * t25 + t69 * t6 + t71 * t5 + t70 * t14) * g(1);
U  = t1;
