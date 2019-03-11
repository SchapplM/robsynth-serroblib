% Calculate potential energy for
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR13_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR13_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:24:02
% EndTime: 2019-03-09 11:24:03
% DurationCPUTime: 0.74s
% Computational Cost: add. (293->108), mult. (550->115), div. (0->0), fcn. (622->12), ass. (0->50)
t71 = -m(1) - m(2);
t70 = -m(5) - m(6);
t69 = mrSges(4,2) - mrSges(3,1);
t68 = mrSges(3,3) + mrSges(4,1);
t67 = -mrSges(4,3) + mrSges(3,2);
t66 = m(6) * qJ(5) - m(7) * (-pkin(10) - qJ(5)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t27 = pkin(11) + qJ(6);
t21 = sin(t27);
t22 = cos(t27);
t28 = sin(pkin(11));
t30 = cos(pkin(11));
t63 = pkin(5) * t28;
t65 = -m(7) * (pkin(9) + t63) - t28 * mrSges(6,1) - t21 * mrSges(7,1) - t30 * mrSges(6,2) - t22 * mrSges(7,2) - mrSges(5,3) + t69;
t20 = pkin(5) * t30 + pkin(4);
t64 = -m(6) * pkin(4) - m(7) * t20 - t30 * mrSges(6,1) - t22 * mrSges(7,1) + t28 * mrSges(6,2) + t21 * mrSges(7,2) - mrSges(5,1);
t29 = sin(pkin(6));
t34 = sin(qJ(2));
t62 = t29 * t34;
t35 = sin(qJ(1));
t61 = t29 * t35;
t37 = cos(qJ(2));
t60 = t29 * t37;
t38 = cos(qJ(1));
t59 = t29 * t38;
t58 = t34 * t35;
t57 = t34 * t38;
t56 = t35 * t37;
t55 = t37 * t38;
t54 = qJ(3) * t37;
t53 = pkin(7) + r_base(3);
t52 = pkin(8) * t59;
t51 = t35 * pkin(1) + r_base(2);
t31 = cos(pkin(6));
t49 = t31 * pkin(8) + t53;
t48 = t38 * pkin(1) + pkin(8) * t61 + r_base(1);
t47 = pkin(2) * t62 + t49;
t10 = -t31 * t55 + t58;
t11 = t31 * t57 + t56;
t46 = t11 * pkin(2) + t10 * qJ(3) + t51;
t45 = t31 * pkin(3) + pkin(9) * t62 + t47;
t12 = t31 * t56 + t57;
t13 = -t31 * t58 + t55;
t44 = t13 * pkin(2) + qJ(3) * t12 + t48;
t43 = pkin(3) * t61 + t44;
t42 = -t29 * t54 + t45;
t40 = (-pkin(3) - pkin(8)) * t59 + t46;
t36 = cos(qJ(4));
t33 = sin(qJ(4));
t9 = t31 * t36 - t33 * t60;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t53 - mrSges(2,3) - m(3) * t49 - m(4) * t47 - m(5) * t42 - t9 * mrSges(5,1) - mrSges(5,3) * t62 - m(6) * (pkin(4) * t9 + t42) - (t28 * t62 + t30 * t9) * mrSges(6,1) - (-t28 * t9 + t30 * t62) * mrSges(6,2) - m(7) * (t20 * t9 + t45) - (t21 * t62 + t22 * t9) * mrSges(7,1) - (-t21 * t9 + t22 * t62) * mrSges(7,2) - t68 * t31 + (m(7) * t54 + (m(4) * qJ(3) - t67) * t37 + (-m(7) * t63 + t69) * t34) * t29 - t66 * (t31 * t33 + t36 * t60)) * g(3) + (-m(7) * t40 - m(3) * (t51 - t52) - m(4) * (t46 - t52) - mrSges(1,2) - t35 * mrSges(2,1) - t38 * mrSges(2,2) + t71 * r_base(2) + t68 * t59 + t70 * (t11 * pkin(9) + t40) + t67 * t10 + t66 * (t10 * t36 + t33 * t59) + t64 * (t10 * t33 - t36 * t59) + t65 * t11) * g(2) + (-m(3) * t48 - m(4) * t44 - m(7) * t43 - t38 * mrSges(2,1) + t35 * mrSges(2,2) - mrSges(1,1) + t71 * r_base(1) - t68 * t61 + t70 * (pkin(9) * t13 + t43) + t67 * t12 + t64 * (t12 * t33 + t36 * t61) + t65 * t13 - t66 * (-t12 * t36 + t33 * t61)) * g(1);
U  = t1;
