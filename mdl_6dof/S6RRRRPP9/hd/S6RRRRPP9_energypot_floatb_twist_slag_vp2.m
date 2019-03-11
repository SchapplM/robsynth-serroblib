% Calculate potential energy for
% S6RRRRPP9
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPP9_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:41
% EndTime: 2019-03-09 21:41:41
% DurationCPUTime: 0.57s
% Computational Cost: add. (356->102), mult. (753->113), div. (0->0), fcn. (903->10), ass. (0->55)
t74 = -m(1) - m(2);
t73 = -mrSges(4,3) + mrSges(3,2);
t72 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t71 = mrSges(4,2) - mrSges(5,3) - mrSges(6,1) - m(7) * (pkin(5) + pkin(10)) - mrSges(7,1);
t70 = -m(7) * qJ(6) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t39 = sin(qJ(2));
t42 = cos(qJ(2));
t43 = cos(qJ(1));
t40 = sin(qJ(1));
t61 = cos(pkin(6));
t57 = t40 * t61;
t26 = -t39 * t57 + t43 * t42;
t38 = sin(qJ(3));
t36 = sin(pkin(6));
t66 = cos(qJ(3));
t59 = t36 * t66;
t14 = t26 * t38 - t40 * t59;
t69 = pkin(10) * t14;
t56 = t43 * t61;
t24 = t39 * t56 + t40 * t42;
t12 = t24 * t38 + t43 * t59;
t68 = t12 * pkin(10);
t65 = t36 * t39;
t21 = t38 * t65 - t61 * t66;
t67 = t21 * pkin(10);
t64 = t36 * t40;
t63 = t36 * t42;
t62 = t36 * t43;
t60 = pkin(7) + r_base(3);
t58 = t61 * pkin(8) + t60;
t55 = t43 * pkin(1) + pkin(8) * t64 + r_base(1);
t53 = t40 * pkin(1) - pkin(8) * t62 + r_base(2);
t25 = t43 * t39 + t42 * t57;
t52 = t26 * pkin(2) + pkin(9) * t25 + t55;
t51 = pkin(2) * t65 - pkin(9) * t63 + t58;
t15 = t26 * t66 + t38 * t64;
t50 = t15 * pkin(3) + t52;
t22 = t38 * t61 + t39 * t59;
t49 = t22 * pkin(3) + t51;
t23 = t39 * t40 - t42 * t56;
t48 = t24 * pkin(2) + t23 * pkin(9) + t53;
t13 = t24 * t66 - t38 * t62;
t47 = t13 * pkin(3) + t48;
t37 = sin(qJ(4));
t41 = cos(qJ(4));
t5 = t15 * t37 - t25 * t41;
t6 = t15 * t41 + t25 * t37;
t46 = t6 * pkin(4) + qJ(5) * t5 + t50;
t10 = t22 * t37 + t41 * t63;
t11 = t22 * t41 - t37 * t63;
t45 = t11 * pkin(4) + qJ(5) * t10 + t49;
t3 = t13 * t37 - t23 * t41;
t4 = t13 * t41 + t23 * t37;
t44 = t4 * pkin(4) + t3 * qJ(5) + t47;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t60 - mrSges(2,3) - m(3) * t58 - t61 * mrSges(3,3) - (t39 * mrSges(3,1) + t42 * mrSges(3,2)) * t36 - m(4) * t51 - t22 * mrSges(4,1) + mrSges(4,3) * t63 - m(5) * (t49 + t67) - m(6) * (t45 + t67) - m(7) * t45 + t70 * t11 + t72 * t10 + t71 * t21) * g(3) + (-m(5) * (t47 + t68) - m(6) * (t44 + t68) + mrSges(3,3) * t62 - m(3) * t53 - m(4) * t48 - m(7) * t44 - mrSges(1,2) - t13 * mrSges(4,1) - t24 * mrSges(3,1) - t40 * mrSges(2,1) - t43 * mrSges(2,2) + t74 * r_base(2) + t73 * t23 + t70 * t4 + t72 * t3 + t71 * t12) * g(2) + (-m(5) * (t50 + t69) - m(6) * (t46 + t69) - m(3) * t55 - m(4) * t52 - m(7) * t46 - mrSges(3,3) * t64 - mrSges(1,1) - t15 * mrSges(4,1) - t26 * mrSges(3,1) + t40 * mrSges(2,2) - t43 * mrSges(2,1) + t74 * r_base(1) + t73 * t25 + t70 * t6 + t72 * t5 + t71 * t14) * g(1);
U  = t1;
