% Calculate potential energy for
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPPR2_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPPR2_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:03:51
% EndTime: 2019-03-08 21:03:52
% DurationCPUTime: 0.72s
% Computational Cost: add. (353->104), mult. (537->115), div. (0->0), fcn. (603->12), ass. (0->51)
t72 = -m(1) - m(2);
t71 = -m(6) - m(7);
t70 = m(4) * pkin(8) - mrSges(3,2) + mrSges(4,3);
t69 = -m(7) * pkin(9) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t39 = sin(qJ(6));
t42 = cos(qJ(6));
t68 = -t39 * mrSges(7,1) - t42 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t67 = m(7) * pkin(5) + t42 * mrSges(7,1) - t39 * mrSges(7,2) + mrSges(6,1) + mrSges(5,3);
t66 = -t67 - t70;
t40 = sin(qJ(3));
t65 = pkin(3) * t40;
t34 = sin(pkin(10));
t35 = sin(pkin(6));
t64 = t34 * t35;
t36 = cos(pkin(10));
t63 = t35 * t36;
t62 = t35 * t40;
t41 = sin(qJ(2));
t61 = t35 * t41;
t43 = cos(qJ(3));
t60 = t35 * t43;
t44 = cos(qJ(2));
t59 = t35 * t44;
t37 = cos(pkin(6));
t58 = t37 * t41;
t57 = t37 * t44;
t56 = t34 * pkin(1) + r_base(2);
t55 = t34 * t62;
t54 = qJ(1) + r_base(3);
t53 = t36 * pkin(1) + pkin(7) * t64 + r_base(1);
t52 = t37 * pkin(7) + t54;
t51 = -pkin(7) * t63 + t56;
t16 = t34 * t57 + t36 * t41;
t17 = -t34 * t58 + t36 * t44;
t27 = pkin(3) * t43 + pkin(2);
t38 = -qJ(4) - pkin(8);
t50 = pkin(3) * t55 - t16 * t38 + t17 * t27 + t53;
t49 = t27 * t61 + t37 * t65 + t38 * t59 + t52;
t14 = t34 * t41 - t36 * t57;
t15 = t34 * t44 + t36 * t58;
t46 = -t14 * t38 + t15 * t27 + (-pkin(7) - t65) * t63 + t56;
t33 = qJ(3) + pkin(11);
t29 = cos(t33);
t28 = sin(t33);
t11 = t28 * t37 + t29 * t61;
t10 = t28 * t61 - t37 * t29;
t6 = t17 * t29 + t28 * t64;
t5 = t17 * t28 - t29 * t64;
t4 = t15 * t29 - t28 * t63;
t3 = t15 * t28 + t29 * t63;
t1 = (-m(1) * r_base(3) - m(2) * t54 - m(5) * t49 - mrSges(1,3) - mrSges(2,3) + (-m(3) - m(4)) * t52 + t71 * (t11 * pkin(4) + t10 * qJ(5) + t49) + (-t40 * mrSges(4,1) - t43 * mrSges(4,2) - mrSges(3,3)) * t37 + (t70 * t44 + (-m(4) * pkin(2) - mrSges(4,1) * t43 + mrSges(4,2) * t40 - mrSges(3,1)) * t41) * t35 + t67 * t59 + t68 * t10 + t69 * t11) * g(3) + (-(t15 * t43 - t36 * t62) * mrSges(4,1) + mrSges(3,3) * t63 - (-t15 * t40 - t36 * t60) * mrSges(4,2) - m(4) * (pkin(2) * t15 + t51) - m(3) * t51 - m(5) * t46 - mrSges(1,2) - t15 * mrSges(3,1) - t34 * mrSges(2,1) - t36 * mrSges(2,2) + t72 * r_base(2) + t71 * (t4 * pkin(4) + qJ(5) * t3 + t46) + t69 * t4 + t68 * t3 + t66 * t14) * g(2) + (-(-t17 * t40 + t34 * t60) * mrSges(4,2) - (t17 * t43 + t55) * mrSges(4,1) - m(4) * (pkin(2) * t17 + t53) - m(3) * t53 - m(5) * t50 - mrSges(3,3) * t64 - mrSges(1,1) - t17 * mrSges(3,1) + t34 * mrSges(2,2) - t36 * mrSges(2,1) + t72 * r_base(1) + t71 * (t6 * pkin(4) + qJ(5) * t5 + t50) + t69 * t6 + t68 * t5 + t66 * t16) * g(1);
U  = t1;
