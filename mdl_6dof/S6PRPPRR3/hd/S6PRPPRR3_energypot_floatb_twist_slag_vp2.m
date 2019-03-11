% Calculate potential energy for
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPPRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:20
% EndTime: 2019-03-08 19:21:21
% DurationCPUTime: 0.71s
% Computational Cost: add. (333->98), mult. (691->118), div. (0->0), fcn. (818->12), ass. (0->49)
t75 = -m(1) - m(2);
t74 = -m(6) - m(7);
t73 = -mrSges(3,1) - mrSges(4,1);
t72 = -mrSges(4,3) + mrSges(3,2);
t71 = mrSges(3,3) + mrSges(4,2) - mrSges(5,3);
t70 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t41 = sin(qJ(6));
t44 = cos(qJ(6));
t69 = -t41 * mrSges(7,1) - t44 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t68 = -m(7) * pkin(5) - t44 * mrSges(7,1) + t41 * mrSges(7,2) - mrSges(6,1);
t36 = sin(pkin(10));
t37 = sin(pkin(6));
t67 = t36 * t37;
t39 = cos(pkin(10));
t66 = t37 * t39;
t42 = sin(qJ(5));
t65 = t37 * t42;
t43 = sin(qJ(2));
t64 = t37 * t43;
t45 = cos(qJ(5));
t63 = t37 * t45;
t40 = cos(pkin(6));
t62 = t40 * t43;
t46 = cos(qJ(2));
t61 = t40 * t46;
t60 = qJ(4) * t37;
t59 = qJ(1) + r_base(3);
t58 = t39 * pkin(1) + pkin(7) * t67 + r_base(1);
t57 = t40 * pkin(7) + t59;
t56 = pkin(2) * t64 + t57;
t55 = t36 * pkin(1) - pkin(7) * t66 + r_base(2);
t24 = t36 * t61 + t39 * t43;
t25 = -t36 * t62 + t39 * t46;
t54 = t25 * pkin(2) + qJ(3) * t24 + t58;
t22 = t36 * t43 - t39 * t61;
t23 = t36 * t46 + t39 * t62;
t53 = t23 * pkin(2) + qJ(3) * t22 + t55;
t52 = t23 * pkin(3) + t39 * t60 + t53;
t51 = t25 * pkin(3) - t36 * t60 + t54;
t50 = -qJ(3) * t37 * t46 + pkin(3) * t64 - t40 * qJ(4) + t56;
t38 = cos(pkin(11));
t35 = sin(pkin(11));
t17 = (-t35 * t46 + t38 * t43) * t37;
t16 = (t35 * t43 + t38 * t46) * t37;
t10 = t24 * t35 + t25 * t38;
t9 = -t24 * t38 + t25 * t35;
t8 = t22 * t35 + t23 * t38;
t7 = -t22 * t38 + t23 * t35;
t1 = (-m(1) * r_base(3) - m(2) * t59 - m(3) * t57 - m(4) * t56 - m(5) * t50 - t17 * mrSges(5,1) - mrSges(1,3) - mrSges(2,3) + t74 * (t17 * pkin(4) + t16 * pkin(8) + t50) + ((m(4) * qJ(3) - t72) * t46 + t73 * t43) * t37 + t68 * (t17 * t45 - t40 * t42) + t69 * t16 + t70 * (t17 * t42 + t40 * t45) - t71 * t40) * g(3) + (-m(3) * t55 - m(4) * t53 - m(5) * t52 - t36 * mrSges(2,1) - t8 * mrSges(5,1) - t39 * mrSges(2,2) - mrSges(1,2) + t75 * r_base(2) + t74 * (t8 * pkin(4) + pkin(8) * t7 + t52) + t68 * (t39 * t65 + t45 * t8) + t69 * t7 + t73 * t23 + t72 * t22 + t70 * (-t39 * t63 + t42 * t8) + t71 * t66) * g(2) + (-m(3) * t58 - m(4) * t54 - m(5) * t51 - t39 * mrSges(2,1) - t10 * mrSges(5,1) + t36 * mrSges(2,2) - mrSges(1,1) + t75 * r_base(1) + t74 * (t10 * pkin(4) + pkin(8) * t9 + t51) + t68 * (t10 * t45 - t36 * t65) + t69 * t9 + t70 * (t10 * t42 + t36 * t63) + t73 * t25 + t72 * t24 - t71 * t67) * g(1);
U  = t1;
