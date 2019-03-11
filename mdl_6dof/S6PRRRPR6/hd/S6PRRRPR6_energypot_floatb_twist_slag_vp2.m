% Calculate potential energy for
% S6PRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRPR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:02
% EndTime: 2019-03-08 23:31:03
% DurationCPUTime: 0.65s
% Computational Cost: add. (382->105), mult. (821->118), div. (0->0), fcn. (997->12), ass. (0->58)
t76 = -m(2) - m(1);
t75 = -mrSges(4,3) + mrSges(3,2);
t74 = mrSges(4,2) - mrSges(5,3) - mrSges(6,2) - m(7) * (pkin(9) - pkin(10)) + mrSges(7,3);
t38 = sin(qJ(6));
t42 = cos(qJ(6));
t73 = -t38 * mrSges(7,1) - t42 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t72 = -m(7) * pkin(5) - t42 * mrSges(7,1) + t38 * mrSges(7,2) - mrSges(5,1) - mrSges(6,1);
t35 = sin(pkin(11));
t37 = cos(pkin(11));
t44 = cos(qJ(2));
t41 = sin(qJ(2));
t62 = cos(pkin(6));
t59 = t41 * t62;
t22 = t35 * t44 + t37 * t59;
t40 = sin(qJ(3));
t36 = sin(pkin(6));
t68 = cos(qJ(3));
t60 = t36 * t68;
t10 = t22 * t40 + t37 * t60;
t71 = pkin(9) * t10;
t24 = -t35 * t59 + t37 * t44;
t12 = t24 * t40 - t35 * t60;
t70 = pkin(9) * t12;
t64 = t36 * t41;
t25 = t40 * t64 - t62 * t68;
t69 = t25 * pkin(9);
t67 = t35 * t36;
t66 = t36 * t37;
t65 = t36 * t40;
t63 = t36 * t44;
t61 = qJ(1) + r_base(3);
t58 = t44 * t62;
t57 = t37 * pkin(1) + pkin(7) * t67 + r_base(1);
t56 = t62 * pkin(7) + t61;
t54 = t35 * pkin(1) - pkin(7) * t66 + r_base(2);
t23 = t35 * t58 + t37 * t41;
t53 = t24 * pkin(2) + pkin(8) * t23 + t57;
t13 = t24 * t68 + t35 * t65;
t52 = t13 * pkin(3) + t53;
t51 = pkin(2) * t64 - pkin(8) * t63 + t56;
t26 = t40 * t62 + t41 * t60;
t50 = t26 * pkin(3) + t51;
t21 = t35 * t41 - t37 * t58;
t49 = t22 * pkin(2) + pkin(8) * t21 + t54;
t11 = t22 * t68 - t37 * t65;
t48 = t11 * pkin(3) + t49;
t39 = sin(qJ(4));
t43 = cos(qJ(4));
t5 = t13 * t39 - t23 * t43;
t6 = t13 * t43 + t23 * t39;
t47 = t6 * pkin(4) + qJ(5) * t5 + t52;
t14 = t26 * t39 + t43 * t63;
t15 = t26 * t43 - t39 * t63;
t46 = t15 * pkin(4) + t14 * qJ(5) + t50;
t3 = t11 * t39 - t21 * t43;
t4 = t11 * t43 + t21 * t39;
t45 = t4 * pkin(4) + qJ(5) * t3 + t48;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t61 - mrSges(2,3) - m(3) * t56 - t62 * mrSges(3,3) - (mrSges(3,1) * t41 + mrSges(3,2) * t44) * t36 - m(4) * t51 - t26 * mrSges(4,1) + mrSges(4,3) * t63 - m(5) * (t50 + t69) - m(6) * (t46 + t69) - m(7) * t46 + t72 * t15 + t73 * t14 + t74 * t25) * g(3) + (-m(7) * t45 - m(3) * t54 + mrSges(3,3) * t66 - m(4) * t49 - m(5) * (t48 + t71) - m(6) * (t45 + t71) - mrSges(1,2) - t11 * mrSges(4,1) - t22 * mrSges(3,1) - t35 * mrSges(2,1) - t37 * mrSges(2,2) + t76 * r_base(2) + t72 * t4 + t73 * t3 + t75 * t21 + t74 * t10) * g(2) + (-m(7) * t47 - m(4) * t53 - m(3) * t57 - mrSges(3,3) * t67 - m(5) * (t52 + t70) - m(6) * (t47 + t70) - mrSges(1,1) - t13 * mrSges(4,1) - t24 * mrSges(3,1) + t35 * mrSges(2,2) - t37 * mrSges(2,1) + t76 * r_base(1) + t72 * t6 + t73 * t5 + t75 * t23 + t74 * t12) * g(1);
U  = t1;
