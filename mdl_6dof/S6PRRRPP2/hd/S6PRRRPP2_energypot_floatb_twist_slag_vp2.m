% Calculate potential energy for
% S6PRRRPP2
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:38
% EndTime: 2019-03-08 22:48:39
% DurationCPUTime: 0.60s
% Computational Cost: add. (356->102), mult. (753->114), div. (0->0), fcn. (903->10), ass. (0->56)
t74 = -m(1) - m(2);
t73 = -mrSges(4,3) + mrSges(3,2);
t72 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t71 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2) - m(7) * (pkin(9) - qJ(6)) + mrSges(7,3);
t70 = -m(7) * pkin(5) - mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t35 = sin(pkin(10));
t37 = cos(pkin(10));
t42 = cos(qJ(2));
t40 = sin(qJ(2));
t60 = cos(pkin(6));
t57 = t40 * t60;
t22 = t35 * t42 + t37 * t57;
t39 = sin(qJ(3));
t36 = sin(pkin(6));
t66 = cos(qJ(3));
t58 = t36 * t66;
t10 = t22 * t39 + t37 * t58;
t69 = pkin(9) * t10;
t24 = -t35 * t57 + t37 * t42;
t12 = t24 * t39 - t35 * t58;
t68 = pkin(9) * t12;
t62 = t36 * t40;
t25 = t39 * t62 - t60 * t66;
t67 = t25 * pkin(9);
t65 = t35 * t36;
t64 = t36 * t37;
t63 = t36 * t39;
t61 = t36 * t42;
t59 = qJ(1) + r_base(3);
t56 = t42 * t60;
t55 = t37 * pkin(1) + pkin(7) * t65 + r_base(1);
t54 = t60 * pkin(7) + t59;
t52 = t35 * pkin(1) - pkin(7) * t64 + r_base(2);
t23 = t35 * t56 + t37 * t40;
t51 = t24 * pkin(2) + pkin(8) * t23 + t55;
t13 = t24 * t66 + t35 * t63;
t50 = t13 * pkin(3) + t51;
t49 = pkin(2) * t62 - pkin(8) * t61 + t54;
t26 = t39 * t60 + t40 * t58;
t48 = t26 * pkin(3) + t49;
t21 = t35 * t40 - t37 * t56;
t47 = t22 * pkin(2) + pkin(8) * t21 + t52;
t11 = t22 * t66 - t37 * t63;
t46 = t11 * pkin(3) + t47;
t38 = sin(qJ(4));
t41 = cos(qJ(4));
t5 = t13 * t38 - t23 * t41;
t6 = t13 * t41 + t23 * t38;
t45 = t6 * pkin(4) + qJ(5) * t5 + t50;
t14 = t26 * t38 + t41 * t61;
t15 = t26 * t41 - t38 * t61;
t44 = t15 * pkin(4) + t14 * qJ(5) + t48;
t3 = t11 * t38 - t21 * t41;
t4 = t11 * t41 + t21 * t38;
t43 = t4 * pkin(4) + qJ(5) * t3 + t46;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t59 - mrSges(2,3) - m(3) * t54 - t60 * mrSges(3,3) - (t40 * mrSges(3,1) + t42 * mrSges(3,2)) * t36 - m(4) * t49 - t26 * mrSges(4,1) + mrSges(4,3) * t61 - m(5) * (t48 + t67) - m(6) * (t44 + t67) - m(7) * t44 + t70 * t15 + t72 * t14 + t71 * t25) * g(3) + (-m(5) * (t46 + t69) - m(6) * (t43 + t69) + mrSges(3,3) * t64 - m(3) * t52 - m(4) * t47 - m(7) * t43 - mrSges(1,2) - t11 * mrSges(4,1) - t22 * mrSges(3,1) - t35 * mrSges(2,1) - t37 * mrSges(2,2) + t74 * r_base(2) + t73 * t21 + t70 * t4 + t72 * t3 + t71 * t10) * g(2) + (-m(5) * (t50 + t68) - m(6) * (t45 + t68) - m(3) * t55 - m(4) * t51 - m(7) * t45 - mrSges(1,1) - mrSges(3,3) * t65 - t13 * mrSges(4,1) - t24 * mrSges(3,1) + t35 * mrSges(2,2) - t37 * mrSges(2,1) + t74 * r_base(1) + t73 * t23 + t70 * t6 + t72 * t5 + t71 * t12) * g(1);
U  = t1;
