% Calculate potential energy for
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:04
% EndTime: 2019-03-08 20:45:05
% DurationCPUTime: 0.76s
% Computational Cost: add. (293->98), mult. (550->106), div. (0->0), fcn. (622->12), ass. (0->48)
t72 = -m(1) - m(2);
t71 = -m(5) - m(6);
t70 = -mrSges(3,1) + mrSges(4,2);
t69 = -mrSges(3,3) - mrSges(4,1);
t68 = -mrSges(4,3) + mrSges(3,2);
t67 = m(6) * pkin(9) - m(7) * (-pkin(10) - pkin(9)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t27 = qJ(5) + qJ(6);
t21 = sin(t27);
t22 = cos(t27);
t35 = cos(qJ(5));
t66 = -t21 * mrSges(7,1) - t35 * mrSges(6,2) - t22 * mrSges(7,2) - mrSges(5,3);
t32 = sin(qJ(5));
t65 = -m(7) * (pkin(5) * t32 + pkin(8)) - t32 * mrSges(6,1) + t66 + t70;
t64 = -m(6) * pkin(4) - m(7) * (pkin(5) * t35 + pkin(4)) - mrSges(6,1) * t35 - mrSges(7,1) * t22 + mrSges(6,2) * t32 + mrSges(7,2) * t21 - mrSges(5,1);
t28 = sin(pkin(11));
t29 = sin(pkin(6));
t63 = t28 * t29;
t30 = cos(pkin(11));
t62 = t29 * t30;
t33 = sin(qJ(4));
t61 = t29 * t33;
t34 = sin(qJ(2));
t60 = t29 * t34;
t36 = cos(qJ(4));
t59 = t29 * t36;
t37 = cos(qJ(2));
t58 = t29 * t37;
t31 = cos(pkin(6));
t57 = t31 * t34;
t56 = t31 * t37;
t55 = t32 * t34;
t54 = qJ(3) * t37;
t53 = pkin(7) * t62;
t52 = t28 * pkin(1) + r_base(2);
t51 = qJ(1) + r_base(3);
t49 = t30 * pkin(1) + pkin(7) * t63 + r_base(1);
t48 = t31 * pkin(7) + t51;
t47 = pkin(2) * t60 + t48;
t8 = t28 * t34 - t30 * t56;
t9 = t28 * t37 + t30 * t57;
t46 = t9 * pkin(2) + t8 * qJ(3) + t52;
t45 = t31 * pkin(3) + pkin(8) * t60 + t47;
t10 = t28 * t56 + t30 * t34;
t11 = -t28 * t57 + t30 * t37;
t44 = t11 * pkin(2) + t10 * qJ(3) + t49;
t43 = pkin(3) * t63 + t44;
t40 = (-pkin(3) - pkin(7)) * t62 + t46;
t1 = (-m(1) * r_base(3) - m(2) * t51 - m(3) * t48 - m(4) * t47 - m(7) * t45 - mrSges(1,3) - mrSges(2,3) + t71 * (-t54 * t29 + t45) + t69 * t31 + t66 * t60 + t64 * (t31 * t36 - t33 * t58) + (-t55 * mrSges(6,1) - m(7) * (pkin(5) * t55 - t54) + (m(4) * qJ(3) - t68) * t37 + t70 * t34) * t29 - t67 * (t31 * t33 + t36 * t58)) * g(3) + (-m(3) * (t52 - t53) - m(4) * (t46 - t53) - m(7) * t40 - mrSges(1,2) - t28 * mrSges(2,1) - t30 * mrSges(2,2) + t72 * r_base(2) + t68 * t8 - t69 * t62 + t71 * (t9 * pkin(8) + t40) + t64 * (-t30 * t59 + t33 * t8) + t65 * t9 + t67 * (t30 * t61 + t36 * t8)) * g(2) + (-m(3) * t49 - m(4) * t44 - m(7) * t43 - t30 * mrSges(2,1) + t28 * mrSges(2,2) - mrSges(1,1) + t72 * r_base(1) + t69 * t63 + t71 * (t11 * pkin(8) + t43) + t68 * t10 + t64 * (t10 * t33 + t28 * t59) + t65 * t11 - t67 * (-t10 * t36 + t28 * t61)) * g(1);
U  = t1;
