% Calculate potential energy for
% S6PRRPRR7
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:08
% EndTime: 2019-03-08 22:30:09
% DurationCPUTime: 0.70s
% Computational Cost: add. (322->98), mult. (626->103), div. (0->0), fcn. (727->12), ass. (0->54)
t71 = -m(1) - m(2);
t66 = pkin(4) + pkin(8);
t70 = -m(6) * pkin(9) + m(7) * (-pkin(10) - pkin(9)) - mrSges(4,1) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t28 = qJ(5) + qJ(6);
t23 = sin(t28);
t24 = cos(t28);
t32 = sin(qJ(5));
t35 = cos(qJ(5));
t69 = -m(7) * (pkin(5) * t32 + qJ(4)) - t32 * mrSges(6,1) - t23 * mrSges(7,1) - t35 * mrSges(6,2) - t24 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3);
t68 = m(6) * t66 + m(7) * (pkin(5) * t35 + t66) + t35 * mrSges(6,1) + t24 * mrSges(7,1) - t32 * mrSges(6,2) - t23 * mrSges(7,2) + mrSges(5,1) + mrSges(4,3);
t67 = mrSges(3,2) - t68;
t29 = sin(pkin(11));
t31 = cos(pkin(11));
t34 = sin(qJ(2));
t36 = cos(qJ(2));
t56 = cos(pkin(6));
t50 = t36 * t56;
t10 = t29 * t34 - t31 * t50;
t65 = t10 * pkin(8);
t12 = t29 * t50 + t31 * t34;
t64 = t12 * pkin(8);
t62 = cos(qJ(3));
t30 = sin(pkin(6));
t61 = t29 * t30;
t60 = t30 * t31;
t33 = sin(qJ(3));
t59 = t30 * t33;
t58 = t30 * t34;
t57 = t30 * t36;
t55 = pkin(8) * t57;
t54 = qJ(1) + r_base(3);
t53 = t30 * t62;
t51 = t34 * t56;
t49 = t31 * pkin(1) + pkin(7) * t61 + r_base(1);
t48 = t56 * pkin(7) + t54;
t13 = -t29 * t51 + t31 * t36;
t47 = t13 * pkin(2) + t49;
t46 = pkin(2) * t58 + t48;
t6 = t13 * t62 + t29 * t59;
t45 = t6 * pkin(3) + t47;
t15 = t56 * t33 + t34 * t53;
t44 = t15 * pkin(3) + t46;
t43 = t29 * pkin(1) - pkin(7) * t60 + r_base(2);
t11 = t29 * t36 + t31 * t51;
t42 = t11 * pkin(2) + t43;
t4 = t11 * t62 - t31 * t59;
t41 = t4 * pkin(3) + t42;
t5 = t13 * t33 - t29 * t53;
t40 = t5 * qJ(4) + t45;
t14 = t33 * t58 - t56 * t62;
t39 = t14 * qJ(4) + t44;
t3 = t11 * t33 + t31 * t53;
t38 = t3 * qJ(4) + t41;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t54 - mrSges(2,3) - m(3) * t48 - t56 * mrSges(3,3) - (t34 * mrSges(3,1) + t36 * mrSges(3,2)) * t30 - m(4) * (t46 - t55) - m(5) * (t39 - t55) - m(6) * t39 - m(7) * t44 + t68 * t57 + t69 * t14 + t70 * t15) * g(3) + (-m(4) * (t42 + t65) - m(5) * (t38 + t65) - m(6) * t38 - m(7) * t41 + mrSges(3,3) * t60 - m(3) * t43 - mrSges(1,2) - t11 * mrSges(3,1) - t29 * mrSges(2,1) - t31 * mrSges(2,2) + t71 * r_base(2) + t69 * t3 + t67 * t10 + t70 * t4) * g(2) + (-m(5) * (t40 + t64) - m(6) * t40 - m(7) * t45 - m(4) * (t47 + t64) - m(3) * t49 - mrSges(1,1) - mrSges(3,3) * t61 - t13 * mrSges(3,1) + t29 * mrSges(2,2) - t31 * mrSges(2,1) + t71 * r_base(1) + t69 * t5 + t67 * t12 + t70 * t6) * g(1);
U  = t1;
