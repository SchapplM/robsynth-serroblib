% Calculate potential energy for
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRP2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:07
% EndTime: 2019-03-08 21:28:07
% DurationCPUTime: 0.74s
% Computational Cost: add. (391->111), mult. (605->129), div. (0->0), fcn. (697->12), ass. (0->51)
t79 = -m(1) - m(2);
t78 = -m(6) - m(7);
t77 = mrSges(5,2) - mrSges(6,3) - mrSges(7,2);
t76 = m(4) * pkin(8) - mrSges(3,2) + mrSges(4,3);
t75 = -m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1);
t74 = -m(7) * qJ(6) + mrSges(6,2) - mrSges(7,3);
t73 = -mrSges(5,3) - t76;
t47 = sin(qJ(3));
t72 = pkin(3) * t47;
t41 = sin(pkin(10));
t42 = sin(pkin(6));
t71 = t41 * t42;
t43 = cos(pkin(10));
t70 = t42 * t43;
t69 = t42 * t47;
t48 = sin(qJ(2));
t68 = t42 * t48;
t50 = cos(qJ(3));
t67 = t42 * t50;
t51 = cos(qJ(2));
t66 = t42 * t51;
t44 = cos(pkin(6));
t65 = t44 * t48;
t64 = t44 * t51;
t63 = t41 * pkin(1) + r_base(2);
t62 = t41 * t69;
t61 = qJ(1) + r_base(3);
t60 = t43 * pkin(1) + pkin(7) * t71 + r_base(1);
t59 = t44 * pkin(7) + t61;
t58 = -pkin(7) * t70 + t63;
t24 = t41 * t64 + t43 * t48;
t25 = -t41 * t65 + t43 * t51;
t34 = pkin(3) * t50 + pkin(2);
t45 = -qJ(4) - pkin(8);
t57 = pkin(3) * t62 - t24 * t45 + t25 * t34 + t60;
t56 = t34 * t68 + t44 * t72 + t45 * t66 + t59;
t22 = t41 * t48 - t43 * t64;
t23 = t41 * t51 + t43 * t65;
t53 = t23 * t34 - t22 * t45 + (-pkin(7) - t72) * t70 + t63;
t49 = cos(qJ(5));
t46 = sin(qJ(5));
t40 = qJ(3) + pkin(11);
t36 = cos(t40);
t35 = sin(t40);
t17 = t35 * t44 + t36 * t68;
t16 = t35 * t68 - t44 * t36;
t10 = t25 * t36 + t35 * t71;
t9 = t25 * t35 - t36 * t71;
t8 = t23 * t36 - t35 * t70;
t7 = t23 * t35 + t36 * t70;
t1 = (-m(1) * r_base(3) - m(2) * t61 - m(5) * t56 - t17 * mrSges(5,1) + mrSges(5,3) * t66 - mrSges(1,3) - mrSges(2,3) + (-m(3) - m(4)) * t59 + t78 * (t17 * pkin(4) + pkin(9) * t16 + t56) + (-t47 * mrSges(4,1) - t50 * mrSges(4,2) - mrSges(3,3)) * t44 + (t76 * t51 + (-m(4) * pkin(2) - t50 * mrSges(4,1) + t47 * mrSges(4,2) - mrSges(3,1)) * t48) * t42 + t75 * (t17 * t49 - t46 * t66) + t74 * (t17 * t46 + t49 * t66) + t77 * t16) * g(3) + (-m(5) * t53 - m(4) * (pkin(2) * t23 + t58) - m(3) * t58 + mrSges(3,3) * t70 - mrSges(1,2) - t8 * mrSges(5,1) - t23 * mrSges(3,1) - t41 * mrSges(2,1) - t43 * mrSges(2,2) - (-t23 * t47 - t43 * t67) * mrSges(4,2) - (t23 * t50 - t43 * t69) * mrSges(4,1) + t79 * r_base(2) + t78 * (t8 * pkin(4) + pkin(9) * t7 + t53) + t75 * (t22 * t46 + t49 * t8) + t74 * (-t22 * t49 + t46 * t8) + t77 * t7 + t73 * t22) * g(2) + (-(t25 * t50 + t62) * mrSges(4,1) - m(4) * (pkin(2) * t25 + t60) - m(5) * t57 - m(3) * t60 - mrSges(3,3) * t71 - mrSges(1,1) - t10 * mrSges(5,1) - t25 * mrSges(3,1) + t41 * mrSges(2,2) - t43 * mrSges(2,1) - (-t25 * t47 + t41 * t67) * mrSges(4,2) + t79 * r_base(1) + t78 * (t10 * pkin(4) + pkin(9) * t9 + t57) + t75 * (t10 * t49 + t24 * t46) + t74 * (t10 * t46 - t24 * t49) + t77 * t9 + t73 * t24) * g(1);
U  = t1;
