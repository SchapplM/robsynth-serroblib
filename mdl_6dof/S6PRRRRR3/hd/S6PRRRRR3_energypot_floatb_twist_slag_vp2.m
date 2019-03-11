% Calculate potential energy for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRRRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRRRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:47:45
% EndTime: 2019-03-09 00:47:46
% DurationCPUTime: 0.77s
% Computational Cost: add. (356->112), mult. (618->125), div. (0->0), fcn. (717->14), ass. (0->49)
t69 = -m(1) - m(2);
t68 = -m(4) - m(5);
t67 = -m(6) - m(7);
t35 = sin(qJ(4));
t53 = pkin(4) * t35 + pkin(8);
t40 = -pkin(10) - pkin(9);
t66 = -m(5) * pkin(9) + m(6) * t40 + m(7) * (-pkin(11) + t40) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3) - mrSges(7,3);
t31 = qJ(4) + qJ(5);
t28 = qJ(6) + t31;
t20 = sin(t28);
t21 = cos(t28);
t23 = sin(t31);
t24 = cos(t31);
t38 = cos(qJ(4));
t62 = pkin(5) * t23 + t53;
t65 = -m(6) * t53 - m(7) * t62 - t35 * mrSges(5,1) - t23 * mrSges(6,1) - t20 * mrSges(7,1) - t38 * mrSges(5,2) - t24 * mrSges(6,2) - t21 * mrSges(7,2) + mrSges(3,2) - mrSges(4,3);
t22 = t38 * pkin(4) + pkin(3);
t13 = pkin(5) * t24 + t22;
t64 = -m(5) * pkin(3) - m(6) * t22 - m(7) * t13 - t38 * mrSges(5,1) - t24 * mrSges(6,1) - t21 * mrSges(7,1) + t35 * mrSges(5,2) + t23 * mrSges(6,2) + t20 * mrSges(7,2) - mrSges(4,1);
t61 = cos(qJ(3));
t32 = sin(pkin(12));
t33 = sin(pkin(6));
t60 = t32 * t33;
t34 = cos(pkin(12));
t59 = t33 * t34;
t36 = sin(qJ(3));
t58 = t33 * t36;
t37 = sin(qJ(2));
t57 = t33 * t37;
t39 = cos(qJ(2));
t56 = t33 * t39;
t55 = cos(pkin(6));
t54 = qJ(1) + r_base(3);
t52 = t33 * t61;
t51 = t37 * t55;
t50 = t39 * t55;
t49 = t34 * pkin(1) + pkin(7) * t60 + r_base(1);
t48 = t55 * pkin(7) + t54;
t10 = -t32 * t51 + t34 * t39;
t47 = t10 * pkin(2) + t49;
t46 = pkin(2) * t57 + t48;
t45 = t32 * pkin(1) - pkin(7) * t59 + r_base(2);
t8 = t32 * t39 + t34 * t51;
t43 = t8 * pkin(2) + t45;
t42 = -pkin(8) * t56 + t46;
t12 = t36 * t55 + t37 * t52;
t9 = t32 * t50 + t34 * t37;
t7 = t32 * t37 - t34 * t50;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t54 - mrSges(2,3) - m(3) * t48 - t55 * mrSges(3,3) - (t37 * mrSges(3,1) + t39 * mrSges(3,2)) * t33 - m(4) * t42 - t12 * mrSges(4,1) + mrSges(4,3) * t56 - m(5) * (pkin(3) * t12 + t42) - (t12 * t38 - t35 * t56) * mrSges(5,1) - (-t12 * t35 - t38 * t56) * mrSges(5,2) - m(6) * (t12 * t22 - t53 * t56 + t46) - (t12 * t24 - t23 * t56) * mrSges(6,1) - (-t12 * t23 - t24 * t56) * mrSges(6,2) - m(7) * (t12 * t13 - t56 * t62 + t46) - (t12 * t21 - t20 * t56) * mrSges(7,1) - (-t12 * t20 - t21 * t56) * mrSges(7,2) + t66 * (t36 * t57 - t55 * t61)) * g(3) + (-m(3) * t45 - t32 * mrSges(2,1) - t8 * mrSges(3,1) - t34 * mrSges(2,2) + mrSges(3,3) * t59 - mrSges(1,2) + t69 * r_base(2) + t67 * t43 + t68 * (t7 * pkin(8) + t43) + t64 * (-t34 * t58 + t61 * t8) + t65 * t7 + t66 * (t34 * t52 + t8 * t36)) * g(2) + (-m(3) * t49 - t34 * mrSges(2,1) - t10 * mrSges(3,1) + t32 * mrSges(2,2) - mrSges(3,3) * t60 - mrSges(1,1) + t69 * r_base(1) + t67 * t47 + t68 * (t9 * pkin(8) + t47) + t64 * (t10 * t61 + t32 * t58) + t65 * t9 + t66 * (t10 * t36 - t32 * t52)) * g(1);
U  = t1;
