% Calculate potential energy for
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRRR13_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:44:50
% EndTime: 2019-03-09 14:44:51
% DurationCPUTime: 0.73s
% Computational Cost: add. (293->98), mult. (550->102), div. (0->0), fcn. (622->12), ass. (0->48)
t72 = -m(1) - m(2);
t71 = -m(5) - m(6);
t70 = -mrSges(3,1) + mrSges(4,2);
t69 = mrSges(3,3) + mrSges(4,1);
t68 = -mrSges(4,3) + mrSges(3,2);
t67 = m(6) * pkin(10) - m(7) * (-pkin(11) - pkin(10)) - mrSges(5,2) + mrSges(6,3) + mrSges(7,3);
t27 = qJ(5) + qJ(6);
t21 = sin(t27);
t22 = cos(t27);
t34 = cos(qJ(5));
t66 = -t21 * mrSges(7,1) - t34 * mrSges(6,2) - t22 * mrSges(7,2) - mrSges(5,3);
t30 = sin(qJ(5));
t65 = -m(7) * (pkin(5) * t30 + pkin(9)) - t30 * mrSges(6,1) + t66 + t70;
t64 = -m(6) * pkin(4) - m(7) * (pkin(5) * t34 + pkin(4)) - t34 * mrSges(6,1) - t22 * mrSges(7,1) + t30 * mrSges(6,2) + t21 * mrSges(7,2) - mrSges(5,1);
t28 = sin(pkin(6));
t32 = sin(qJ(2));
t63 = t28 * t32;
t33 = sin(qJ(1));
t62 = t28 * t33;
t36 = cos(qJ(2));
t61 = t28 * t36;
t37 = cos(qJ(1));
t60 = t28 * t37;
t59 = t30 * t32;
t58 = t32 * t33;
t57 = t32 * t37;
t56 = t33 * t36;
t55 = t36 * t37;
t54 = qJ(3) * t36;
t53 = pkin(7) + r_base(3);
t52 = pkin(8) * t60;
t51 = t33 * pkin(1) + r_base(2);
t29 = cos(pkin(6));
t49 = t29 * pkin(8) + t53;
t48 = t37 * pkin(1) + pkin(8) * t62 + r_base(1);
t47 = pkin(2) * t63 + t49;
t10 = -t29 * t55 + t58;
t11 = t29 * t57 + t56;
t46 = t11 * pkin(2) + t10 * qJ(3) + t51;
t45 = t29 * pkin(3) + pkin(9) * t63 + t47;
t12 = t29 * t56 + t57;
t13 = -t29 * t58 + t55;
t44 = t13 * pkin(2) + t12 * qJ(3) + t48;
t43 = pkin(3) * t62 + t44;
t40 = (-pkin(3) - pkin(8)) * t60 + t46;
t35 = cos(qJ(4));
t31 = sin(qJ(4));
t1 = (-m(1) * r_base(3) - m(2) * t53 - m(3) * t49 - m(4) * t47 - m(7) * t45 - mrSges(1,3) - mrSges(2,3) + t71 * (-t28 * t54 + t45) + t66 * t63 + t64 * (t29 * t35 - t31 * t61) - t69 * t29 + (-t59 * mrSges(6,1) - m(7) * (pkin(5) * t59 - t54) + (m(4) * qJ(3) - t68) * t36 + t70 * t32) * t28 - t67 * (t29 * t31 + t35 * t61)) * g(3) + (-m(3) * (t51 - t52) - m(4) * (t46 - t52) - m(7) * t40 - mrSges(1,2) - t33 * mrSges(2,1) - t37 * mrSges(2,2) + t72 * r_base(2) + t69 * t60 + t71 * (t11 * pkin(9) + t40) + t68 * t10 + t67 * (t10 * t35 + t31 * t60) + t64 * (t10 * t31 - t35 * t60) + t65 * t11) * g(2) + (-m(3) * t48 - m(4) * t44 - m(7) * t43 - t37 * mrSges(2,1) + t33 * mrSges(2,2) - mrSges(1,1) + t72 * r_base(1) - t69 * t62 + t71 * (t13 * pkin(9) + t43) + t68 * t12 + t64 * (t12 * t31 + t35 * t62) + t65 * t13 - t67 * (-t12 * t35 + t31 * t62)) * g(1);
U  = t1;
