% Calculate potential energy for
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:42:51
% EndTime: 2019-03-08 19:42:52
% DurationCPUTime: 0.74s
% Computational Cost: add. (353->104), mult. (537->113), div. (0->0), fcn. (603->12), ass. (0->49)
t70 = -m(1) - m(2);
t69 = -m(6) - m(7);
t68 = m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3);
t67 = -m(7) * pkin(9) - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t41 = sin(qJ(6));
t43 = cos(qJ(6));
t66 = -t41 * mrSges(7,1) - t43 * mrSges(7,2) + mrSges(5,2) - mrSges(6,3);
t65 = m(7) * pkin(5) + t43 * mrSges(7,1) - t41 * mrSges(7,2) + mrSges(6,1) + mrSges(5,3);
t64 = -t65 - t68;
t34 = sin(pkin(11));
t63 = pkin(3) * t34;
t35 = sin(pkin(10));
t36 = sin(pkin(6));
t62 = t35 * t36;
t38 = cos(pkin(10));
t61 = t36 * t38;
t42 = sin(qJ(2));
t60 = t36 * t42;
t44 = cos(qJ(2));
t59 = t36 * t44;
t39 = cos(pkin(6));
t58 = t39 * t42;
t57 = t39 * t44;
t56 = t35 * pkin(1) + r_base(2);
t55 = t34 * t62;
t54 = qJ(1) + r_base(3);
t53 = t38 * pkin(1) + pkin(7) * t62 + r_base(1);
t52 = t39 * pkin(7) + t54;
t51 = -pkin(7) * t61 + t56;
t16 = t35 * t57 + t38 * t42;
t17 = -t35 * t58 + t38 * t44;
t37 = cos(pkin(11));
t27 = pkin(3) * t37 + pkin(2);
t40 = -pkin(8) - qJ(3);
t50 = pkin(3) * t55 - t16 * t40 + t17 * t27 + t53;
t49 = t27 * t60 + t39 * t63 + t40 * t59 + t52;
t14 = t35 * t42 - t38 * t57;
t15 = t35 * t44 + t38 * t58;
t46 = -t14 * t40 + t15 * t27 + (-pkin(7) - t63) * t61 + t56;
t33 = pkin(11) + qJ(4);
t29 = cos(t33);
t28 = sin(t33);
t11 = t28 * t39 + t29 * t60;
t10 = t28 * t60 - t39 * t29;
t6 = t17 * t29 + t28 * t62;
t5 = t17 * t28 - t29 * t62;
t4 = t15 * t29 - t28 * t61;
t3 = t15 * t28 + t29 * t61;
t1 = (-m(1) * r_base(3) - m(2) * t54 - m(5) * t49 - mrSges(1,3) - mrSges(2,3) + (-m(3) - m(4)) * t52 + t69 * (t11 * pkin(4) + t10 * qJ(5) + t49) + (-t34 * mrSges(4,1) - t37 * mrSges(4,2) - mrSges(3,3)) * t39 + (t68 * t44 + (-m(4) * pkin(2) - mrSges(4,1) * t37 + mrSges(4,2) * t34 - mrSges(3,1)) * t42) * t36 + t65 * t59 + t66 * t10 + t67 * t11) * g(3) + (-m(3) * t51 + mrSges(3,3) * t61 - m(5) * t46 - mrSges(1,2) - t15 * mrSges(3,1) - t35 * mrSges(2,1) - t38 * mrSges(2,2) - m(4) * (pkin(2) * t15 + t51) - (-t15 * t34 - t37 * t61) * mrSges(4,2) - (t15 * t37 - t34 * t61) * mrSges(4,1) + t70 * r_base(2) + t69 * (t4 * pkin(4) + qJ(5) * t3 + t46) + t67 * t4 + t66 * t3 + t64 * t14) * g(2) + (-m(3) * t53 - mrSges(3,3) * t62 - m(5) * t50 - mrSges(1,1) - t17 * mrSges(3,1) + t35 * mrSges(2,2) - t38 * mrSges(2,1) - m(4) * (pkin(2) * t17 + t53) - (t17 * t37 + t55) * mrSges(4,1) - (-t17 * t34 + t37 * t62) * mrSges(4,2) + t70 * r_base(1) + t69 * (t6 * pkin(4) + qJ(5) * t5 + t50) + t67 * t6 + t66 * t5 + t64 * t16) * g(1);
U  = t1;
