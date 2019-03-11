% Calculate potential energy for
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:24
% EndTime: 2019-03-08 21:12:25
% DurationCPUTime: 0.65s
% Computational Cost: add. (382->105), mult. (821->118), div. (0->0), fcn. (997->12), ass. (0->58)
t76 = -m(2) - m(1);
t75 = -mrSges(4,3) + mrSges(3,2);
t74 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2) - m(7) * (-pkin(9) + qJ(4)) + mrSges(7,3);
t40 = sin(qJ(6));
t43 = cos(qJ(6));
t73 = -t40 * mrSges(7,1) - t43 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t72 = -m(7) * pkin(5) - t43 * mrSges(7,1) + t40 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t71 = cos(qJ(3));
t36 = sin(pkin(10));
t37 = sin(pkin(6));
t70 = t36 * t37;
t39 = cos(pkin(10));
t69 = t37 * t39;
t41 = sin(qJ(3));
t68 = t37 * t41;
t42 = sin(qJ(2));
t67 = t37 * t42;
t44 = cos(qJ(2));
t66 = t37 * t44;
t62 = cos(pkin(6));
t59 = t42 * t62;
t22 = t36 * t44 + t39 * t59;
t60 = t37 * t71;
t12 = t22 * t41 + t39 * t60;
t65 = qJ(4) * t12;
t24 = -t36 * t59 + t39 * t44;
t14 = t24 * t41 - t36 * t60;
t64 = qJ(4) * t14;
t25 = t41 * t67 - t62 * t71;
t63 = t25 * qJ(4);
t61 = qJ(1) + r_base(3);
t58 = t44 * t62;
t57 = t39 * pkin(1) + pkin(7) * t70 + r_base(1);
t56 = t62 * pkin(7) + t61;
t54 = t36 * pkin(1) - pkin(7) * t69 + r_base(2);
t23 = t36 * t58 + t39 * t42;
t53 = t24 * pkin(2) + pkin(8) * t23 + t57;
t15 = t24 * t71 + t36 * t68;
t52 = t15 * pkin(3) + t53;
t51 = pkin(2) * t67 - pkin(8) * t66 + t56;
t26 = t62 * t41 + t42 * t60;
t50 = t26 * pkin(3) + t51;
t21 = t36 * t42 - t39 * t58;
t49 = t22 * pkin(2) + pkin(8) * t21 + t54;
t13 = t22 * t71 - t39 * t68;
t48 = t13 * pkin(3) + t49;
t35 = sin(pkin(11));
t38 = cos(pkin(11));
t5 = t15 * t35 - t23 * t38;
t6 = t15 * t38 + t23 * t35;
t47 = t6 * pkin(4) + qJ(5) * t5 + t52;
t10 = t26 * t35 + t38 * t66;
t11 = t26 * t38 - t35 * t66;
t46 = t11 * pkin(4) + t10 * qJ(5) + t50;
t3 = t13 * t35 - t21 * t38;
t4 = t13 * t38 + t21 * t35;
t45 = t4 * pkin(4) + qJ(5) * t3 + t48;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t61 - mrSges(2,3) - m(3) * t56 - t62 * mrSges(3,3) - (mrSges(3,1) * t42 + mrSges(3,2) * t44) * t37 - m(4) * t51 - t26 * mrSges(4,1) + mrSges(4,3) * t66 - m(5) * (t50 + t63) - m(6) * (t46 + t63) - m(7) * t46 + t72 * t11 + t73 * t10 + t74 * t25) * g(3) + (mrSges(3,3) * t69 - m(5) * (t48 + t65) - m(6) * (t45 + t65) - m(4) * t49 - m(3) * t54 - m(7) * t45 - mrSges(1,2) - t13 * mrSges(4,1) - t22 * mrSges(3,1) - t36 * mrSges(2,1) - t39 * mrSges(2,2) + t76 * r_base(2) + t72 * t4 + t73 * t3 + t75 * t21 + t74 * t12) * g(2) + (-m(5) * (t52 + t64) - m(6) * (t47 + t64) - m(3) * t57 - m(4) * t53 - m(7) * t47 - mrSges(3,3) * t70 - mrSges(1,1) - t15 * mrSges(4,1) - t24 * mrSges(3,1) + t36 * mrSges(2,2) - t39 * mrSges(2,1) + t76 * r_base(1) + t72 * t6 + t73 * t5 + t75 * t23 + t74 * t14) * g(1);
U  = t1;
