% Calculate potential energy for
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR4_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:09
% EndTime: 2019-03-08 22:10:10
% DurationCPUTime: 0.74s
% Computational Cost: add. (355->101), mult. (751->114), div. (0->0), fcn. (900->12), ass. (0->51)
t80 = -m(6) - m(7);
t40 = sin(qJ(6));
t44 = cos(qJ(6));
t82 = -t40 * mrSges(7,1) - t44 * mrSges(7,2) + mrSges(5,2) + mrSges(4,3) - mrSges(6,3);
t81 = -m(1) - m(2);
t79 = -mrSges(4,1) - mrSges(5,1);
t78 = mrSges(4,2) - mrSges(5,3);
t77 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t75 = -m(7) * pkin(5) - t44 * mrSges(7,1) + t40 * mrSges(7,2) - mrSges(6,1);
t74 = mrSges(3,2) - t82 + t80 * (pkin(8) - pkin(9));
t37 = sin(pkin(11));
t39 = cos(pkin(11));
t43 = sin(qJ(2));
t46 = cos(qJ(2));
t64 = cos(pkin(6));
t60 = t46 * t64;
t22 = t37 * t43 - t39 * t60;
t72 = pkin(8) * t22;
t24 = t37 * t60 + t39 * t43;
t71 = pkin(8) * t24;
t70 = cos(qJ(3));
t38 = sin(pkin(6));
t69 = t37 * t38;
t68 = t38 * t39;
t42 = sin(qJ(3));
t67 = t38 * t42;
t66 = t38 * t43;
t65 = t38 * t46;
t63 = qJ(1) + r_base(3);
t62 = t38 * t70;
t61 = t43 * t64;
t59 = t39 * pkin(1) + pkin(7) * t69 + r_base(1);
t58 = t64 * pkin(7) + t63;
t25 = -t37 * t61 + t39 * t46;
t56 = t25 * pkin(2) + t59;
t55 = t37 * pkin(1) - pkin(7) * t68 + r_base(2);
t23 = t37 * t46 + t39 * t61;
t54 = t23 * pkin(2) + t55;
t53 = pkin(2) * t66 - pkin(8) * t65 + t58;
t15 = t25 * t42 - t37 * t62;
t16 = t25 * t70 + t37 * t67;
t52 = t16 * pkin(3) + qJ(4) * t15 + t56;
t13 = t23 * t42 + t39 * t62;
t14 = t23 * t70 - t39 * t67;
t50 = t14 * pkin(3) + qJ(4) * t13 + t54;
t26 = t42 * t66 - t64 * t70;
t27 = t64 * t42 + t43 * t62;
t48 = t27 * pkin(3) + t26 * qJ(4) + t53;
t45 = cos(qJ(5));
t41 = sin(qJ(5));
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t63 - mrSges(2,3) - m(3) * t58 - t64 * mrSges(3,3) - (t43 * mrSges(3,1) + t46 * mrSges(3,2)) * t38 - m(4) * t53 - m(5) * t48 + t80 * (t27 * pkin(4) + pkin(9) * t65 + t48) + t77 * (-t26 * t45 + t27 * t41) + t79 * t27 + t78 * t26 + t75 * (t26 * t41 + t27 * t45) + t82 * t65) * g(3) + (-mrSges(1,2) - t23 * mrSges(3,1) - t37 * mrSges(2,1) - t39 * mrSges(2,2) - m(3) * t55 + mrSges(3,3) * t68 - m(4) * (t54 + t72) - m(5) * (t50 + t72) + t81 * r_base(2) + t80 * (t14 * pkin(4) + t50) + t79 * t14 + t78 * t13 + t77 * (-t13 * t45 + t14 * t41) + t75 * (t13 * t41 + t14 * t45) + t74 * t22) * g(2) + (-mrSges(1,1) - mrSges(3,3) * t69 - t25 * mrSges(3,1) + t37 * mrSges(2,2) - t39 * mrSges(2,1) - m(3) * t59 - m(4) * (t56 + t71) - m(5) * (t52 + t71) + t81 * r_base(1) + t80 * (t16 * pkin(4) + t52) + t77 * (-t15 * t45 + t16 * t41) + t79 * t16 + t78 * t15 + t75 * (t15 * t41 + t16 * t45) + t74 * t24) * g(1);
U  = t1;
