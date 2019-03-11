% Calculate potential energy for
% S6RRRRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRRPR13_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRRPR13_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR13_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR13_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:51:02
% EndTime: 2019-03-09 23:51:03
% DurationCPUTime: 0.64s
% Computational Cost: add. (382->105), mult. (821->117), div. (0->0), fcn. (997->12), ass. (0->57)
t75 = -m(2) - m(1);
t74 = -mrSges(4,3) + mrSges(3,2);
t73 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2) - m(7) * (pkin(10) - pkin(11)) + mrSges(7,3);
t36 = sin(qJ(6));
t41 = cos(qJ(6));
t72 = -t36 * mrSges(7,1) - t41 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t71 = -m(7) * pkin(5) - t41 * mrSges(7,1) + t36 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t39 = sin(qJ(2));
t43 = cos(qJ(2));
t44 = cos(qJ(1));
t40 = sin(qJ(1));
t62 = cos(pkin(6));
t58 = t40 * t62;
t26 = -t39 * t58 + t44 * t43;
t38 = sin(qJ(3));
t35 = sin(pkin(6));
t67 = cos(qJ(3));
t60 = t35 * t67;
t14 = t26 * t38 - t40 * t60;
t70 = pkin(10) * t14;
t66 = t35 * t39;
t21 = t38 * t66 - t62 * t67;
t69 = pkin(10) * t21;
t57 = t44 * t62;
t24 = t39 * t57 + t40 * t43;
t12 = t24 * t38 + t44 * t60;
t68 = t12 * pkin(10);
t65 = t35 * t40;
t64 = t35 * t43;
t63 = t35 * t44;
t61 = pkin(7) + r_base(3);
t59 = t62 * pkin(8) + t61;
t56 = t44 * pkin(1) + pkin(8) * t65 + r_base(1);
t54 = t40 * pkin(1) - pkin(8) * t63 + r_base(2);
t25 = t44 * t39 + t43 * t58;
t53 = t26 * pkin(2) + pkin(9) * t25 + t56;
t52 = pkin(2) * t66 - pkin(9) * t64 + t59;
t15 = t26 * t67 + t38 * t65;
t51 = t15 * pkin(3) + t53;
t22 = t38 * t62 + t39 * t60;
t50 = t22 * pkin(3) + t52;
t23 = t39 * t40 - t43 * t57;
t49 = t24 * pkin(2) + t23 * pkin(9) + t54;
t13 = t24 * t67 - t38 * t63;
t48 = t13 * pkin(3) + t49;
t37 = sin(qJ(4));
t42 = cos(qJ(4));
t5 = t15 * t37 - t25 * t42;
t6 = t15 * t42 + t25 * t37;
t47 = t6 * pkin(4) + qJ(5) * t5 + t51;
t10 = t22 * t37 + t42 * t64;
t11 = t22 * t42 - t37 * t64;
t46 = t11 * pkin(4) + qJ(5) * t10 + t50;
t3 = t13 * t37 - t23 * t42;
t4 = t13 * t42 + t23 * t37;
t45 = t4 * pkin(4) + t3 * qJ(5) + t48;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t61 - mrSges(2,3) - m(3) * t59 - t62 * mrSges(3,3) - (mrSges(3,1) * t39 + mrSges(3,2) * t43) * t35 - m(4) * t52 - t22 * mrSges(4,1) + mrSges(4,3) * t64 - m(5) * (t50 + t69) - m(6) * (t46 + t69) - m(7) * t46 + t71 * t11 + t72 * t10 + t73 * t21) * g(3) + (-m(5) * (t48 + t68) - m(6) * (t45 + t68) - m(7) * t45 - m(3) * t54 + mrSges(3,3) * t63 - m(4) * t49 - mrSges(1,2) - t13 * mrSges(4,1) - t24 * mrSges(3,1) - t40 * mrSges(2,1) - t44 * mrSges(2,2) + t75 * r_base(2) + t71 * t4 + t72 * t3 + t74 * t23 + t73 * t12) * g(2) + (-m(5) * (t51 + t70) - m(6) * (t47 + t70) - m(7) * t47 - m(3) * t56 - mrSges(3,3) * t65 - m(4) * t53 - mrSges(1,1) - t15 * mrSges(4,1) - t26 * mrSges(3,1) + t40 * mrSges(2,2) - t44 * mrSges(2,1) + t75 * r_base(1) + t71 * t6 + t72 * t5 + t74 * t25 + t73 * t14) * g(1);
U  = t1;
