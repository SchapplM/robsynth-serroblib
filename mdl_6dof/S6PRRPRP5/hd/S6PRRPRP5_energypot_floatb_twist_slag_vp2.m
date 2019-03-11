% Calculate potential energy for
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:31
% EndTime: 2019-03-08 21:45:32
% DurationCPUTime: 0.62s
% Computational Cost: add. (322->97), mult. (664->111), div. (0->0), fcn. (780->10), ass. (0->48)
t74 = -m(1) - m(2);
t73 = -m(6) - m(7);
t72 = mrSges(4,2) - mrSges(5,3);
t71 = mrSges(4,3) + mrSges(5,1);
t70 = mrSges(3,2) - t71;
t69 = -mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t68 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t67 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t66 = cos(qJ(3));
t35 = sin(pkin(10));
t36 = sin(pkin(6));
t65 = t35 * t36;
t37 = cos(pkin(10));
t64 = t36 * t37;
t39 = sin(qJ(3));
t63 = t36 * t39;
t40 = sin(qJ(2));
t62 = t36 * t40;
t42 = cos(qJ(2));
t61 = t36 * t42;
t60 = cos(pkin(6));
t59 = pkin(8) * t61;
t58 = qJ(1) + r_base(3);
t57 = t36 * t66;
t56 = t40 * t60;
t55 = t42 * t60;
t54 = t37 * pkin(1) + pkin(7) * t65 + r_base(1);
t53 = t60 * pkin(7) + t58;
t52 = pkin(2) * t62 + t53;
t51 = t35 * pkin(1) - pkin(7) * t64 + r_base(2);
t22 = t35 * t55 + t37 * t40;
t23 = -t35 * t56 + t37 * t42;
t50 = t23 * pkin(2) + pkin(8) * t22 + t54;
t24 = t39 * t62 - t60 * t66;
t25 = t60 * t39 + t40 * t57;
t49 = t25 * pkin(3) + t24 * qJ(4) + t52;
t20 = t35 * t40 - t37 * t55;
t21 = t35 * t42 + t37 * t56;
t48 = t21 * pkin(2) + pkin(8) * t20 + t51;
t11 = t23 * t39 - t35 * t57;
t12 = t23 * t66 + t35 * t63;
t47 = t12 * pkin(3) + qJ(4) * t11 + t50;
t10 = t21 * t66 - t37 * t63;
t9 = t21 * t39 + t37 * t57;
t46 = t10 * pkin(3) + qJ(4) * t9 + t48;
t41 = cos(qJ(5));
t38 = sin(qJ(5));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t58 - mrSges(2,3) - m(3) * t53 - t60 * mrSges(3,3) - (t40 * mrSges(3,1) + t42 * mrSges(3,2)) * t36 - m(4) * (t52 - t59) - m(5) * (t49 - t59) + t71 * t61 + t73 * (t25 * pkin(9) + (-pkin(4) - pkin(8)) * t61 + t49) + t72 * t24 + t68 * (t24 * t38 - t41 * t61) + t67 * (t24 * t41 + t38 * t61) + t69 * t25) * g(3) + (-m(3) * t51 - m(4) * t48 - m(5) * t46 - t35 * mrSges(2,1) - t21 * mrSges(3,1) - t37 * mrSges(2,2) + mrSges(3,3) * t64 - mrSges(1,2) + t74 * r_base(2) + t72 * t9 + t73 * (t20 * pkin(4) + pkin(9) * t10 + t46) + t68 * (t20 * t41 + t38 * t9) - t67 * (t20 * t38 - t9 * t41) + t70 * t20 + t69 * t10) * g(2) + (-m(3) * t54 - m(4) * t50 - m(5) * t47 - t37 * mrSges(2,1) - t23 * mrSges(3,1) + t35 * mrSges(2,2) - mrSges(3,3) * t65 - mrSges(1,1) + t74 * r_base(1) + t73 * (t22 * pkin(4) + pkin(9) * t12 + t47) + t68 * (t11 * t38 + t22 * t41) - t67 * (-t11 * t41 + t22 * t38) + t72 * t11 + t70 * t22 + t69 * t12) * g(1);
U  = t1;
