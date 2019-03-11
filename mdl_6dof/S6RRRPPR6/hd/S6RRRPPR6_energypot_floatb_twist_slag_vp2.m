% Calculate potential energy for
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR6_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:39
% EndTime: 2019-03-09 15:47:40
% DurationCPUTime: 0.74s
% Computational Cost: add. (353->104), mult. (537->111), div. (0->0), fcn. (603->12), ass. (0->51)
t72 = -m(1) - m(2);
t71 = -m(6) - m(7);
t70 = m(4) * pkin(9) - mrSges(3,2) + mrSges(4,3);
t69 = -m(7) * pkin(10) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t37 = sin(qJ(6));
t41 = cos(qJ(6));
t68 = -t37 * mrSges(7,1) - t41 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t67 = m(7) * pkin(5) + t41 * mrSges(7,1) - t37 * mrSges(7,2) + mrSges(6,1) + mrSges(5,3);
t66 = -t67 - t70;
t38 = sin(qJ(3));
t65 = pkin(3) * t38;
t34 = sin(pkin(6));
t39 = sin(qJ(2));
t64 = t34 * t39;
t40 = sin(qJ(1));
t63 = t34 * t40;
t43 = cos(qJ(2));
t62 = t34 * t43;
t44 = cos(qJ(1));
t61 = t34 * t44;
t60 = t39 * t40;
t59 = t39 * t44;
t58 = t40 * t43;
t57 = t43 * t44;
t56 = pkin(7) + r_base(3);
t55 = t40 * pkin(1) + r_base(2);
t54 = t38 * t63;
t35 = cos(pkin(6));
t53 = t35 * pkin(8) + t56;
t52 = t44 * pkin(1) + pkin(8) * t63 + r_base(1);
t51 = -pkin(8) * t61 + t55;
t42 = cos(qJ(3));
t27 = pkin(3) * t42 + pkin(2);
t36 = -qJ(4) - pkin(9);
t50 = t27 * t64 + t35 * t65 + t36 * t62 + t53;
t16 = t35 * t58 + t59;
t17 = -t35 * t60 + t57;
t49 = pkin(3) * t54 - t16 * t36 + t17 * t27 + t52;
t14 = -t35 * t57 + t60;
t15 = t35 * t59 + t58;
t46 = -t14 * t36 + t15 * t27 + (-pkin(8) - t65) * t61 + t55;
t33 = qJ(3) + pkin(11);
t29 = cos(t33);
t28 = sin(t33);
t11 = t28 * t35 + t29 * t64;
t10 = t28 * t64 - t35 * t29;
t6 = t17 * t29 + t28 * t63;
t5 = t17 * t28 - t29 * t63;
t4 = t15 * t29 - t28 * t61;
t3 = t15 * t28 + t29 * t61;
t1 = (-m(1) * r_base(3) - m(2) * t56 - m(5) * t50 - mrSges(1,3) - mrSges(2,3) + (-m(3) - m(4)) * t53 + t71 * (t11 * pkin(4) + qJ(5) * t10 + t50) + (-t38 * mrSges(4,1) - t42 * mrSges(4,2) - mrSges(3,3)) * t35 + (t70 * t43 + (-m(4) * pkin(2) - mrSges(4,1) * t42 + mrSges(4,2) * t38 - mrSges(3,1)) * t39) * t34 + t67 * t62 + t68 * t10 + t69 * t11) * g(3) + (mrSges(3,3) * t61 - (-t15 * t38 - t42 * t61) * mrSges(4,2) - (t15 * t42 - t38 * t61) * mrSges(4,1) - m(4) * (t15 * pkin(2) + t51) - m(3) * t51 - m(5) * t46 - mrSges(1,2) - t15 * mrSges(3,1) - t40 * mrSges(2,1) - t44 * mrSges(2,2) + t72 * r_base(2) + t71 * (t4 * pkin(4) + t3 * qJ(5) + t46) + t69 * t4 + t68 * t3 + t66 * t14) * g(2) + (-(-t17 * t38 + t42 * t63) * mrSges(4,2) - (t17 * t42 + t54) * mrSges(4,1) - m(4) * (pkin(2) * t17 + t52) - m(3) * t52 - m(5) * t49 - mrSges(3,3) * t63 - mrSges(1,1) - t17 * mrSges(3,1) + t40 * mrSges(2,2) - t44 * mrSges(2,1) + t72 * r_base(1) + t71 * (t6 * pkin(4) + qJ(5) * t5 + t49) + t69 * t6 + t68 * t5 + t66 * t16) * g(1);
U  = t1;
