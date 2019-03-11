% Calculate potential energy for
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRPRPR5_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:28:41
% EndTime: 2019-03-09 10:28:42
% DurationCPUTime: 0.80s
% Computational Cost: add. (421->100), mult. (867->115), div. (0->0), fcn. (1059->14), ass. (0->51)
t77 = -m(5) - m(6);
t76 = -mrSges(3,3) - mrSges(4,3);
t36 = sin(pkin(11));
t42 = sin(qJ(2));
t45 = cos(qJ(2));
t61 = cos(pkin(11));
t21 = -t42 * t36 + t45 * t61;
t75 = -m(3) - m(2) - m(1);
t74 = -m(3) * pkin(1) - mrSges(2,1);
t73 = m(3) * pkin(8) - t76;
t72 = -m(6) * qJ(5) + m(7) * (-pkin(10) - qJ(5)) + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t34 = pkin(12) + qJ(6);
t30 = sin(t34);
t31 = cos(t34);
t35 = sin(pkin(12));
t38 = cos(pkin(12));
t71 = m(7) * (pkin(5) * t35 + pkin(9)) + t35 * mrSges(6,1) + t30 * mrSges(7,1) + t38 * mrSges(6,2) + t31 * mrSges(7,2) - mrSges(4,2) + mrSges(5,3);
t70 = -m(6) * pkin(4) - m(7) * (pkin(5) * t38 + pkin(4)) - t38 * mrSges(6,1) - t31 * mrSges(7,1) + t35 * mrSges(6,2) + t30 * mrSges(7,2) - mrSges(5,1);
t69 = pkin(2) * t42;
t37 = sin(pkin(6));
t43 = sin(qJ(1));
t68 = t37 * t43;
t46 = cos(qJ(1));
t67 = t37 * t46;
t65 = t42 * t46;
t64 = t43 * t42;
t63 = t43 * t45;
t62 = t45 * t46;
t60 = pkin(7) + r_base(3);
t39 = cos(pkin(6));
t58 = t39 * pkin(8) + t60;
t19 = t39 * t69 + (-pkin(8) - qJ(3)) * t37;
t29 = pkin(2) * t45 + pkin(1);
t56 = t46 * t19 + t43 * t29 + r_base(2);
t20 = -t45 * t36 - t42 * t61;
t18 = t20 * t39;
t8 = -t18 * t46 + t43 * t21;
t55 = t8 * pkin(3) + t56;
t54 = -t19 * t43 + t46 * t29 + r_base(1);
t53 = t39 * qJ(3) + t37 * t69 + t58;
t10 = t43 * t18 + t21 * t46;
t52 = t10 * pkin(3) + t54;
t17 = t20 * t37;
t51 = -t17 * pkin(3) + t53;
t48 = t39 * t21;
t44 = cos(qJ(4));
t41 = sin(qJ(4));
t16 = t21 * t37;
t9 = t20 * t46 - t43 * t48;
t7 = t43 * t20 + t46 * t48;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t60 - mrSges(2,3) - m(3) * t58 - (t42 * mrSges(3,1) + t45 * mrSges(3,2)) * t37 - m(4) * t53 + t17 * mrSges(4,1) - m(7) * t51 + t77 * (-pkin(9) * t16 + t51) + t76 * t39 + t70 * (-t17 * t44 + t39 * t41) + t71 * t16 + t72 * (-t17 * t41 - t39 * t44)) * g(3) + (-m(4) * t56 - m(7) * t55 - t46 * mrSges(2,2) - (t39 * t62 - t64) * mrSges(3,2) - (t39 * t65 + t63) * mrSges(3,1) - mrSges(1,2) - t8 * mrSges(4,1) + t74 * t43 + t75 * r_base(2) + t77 * (-pkin(9) * t7 + t55) + t70 * (-t41 * t67 + t8 * t44) + t71 * t7 + t73 * t67 + t72 * (t8 * t41 + t44 * t67)) * g(2) + (-m(7) * t52 - m(4) * t54 + t43 * mrSges(2,2) - (-t39 * t64 + t62) * mrSges(3,1) - (-t39 * t63 - t65) * mrSges(3,2) - mrSges(1,1) - t10 * mrSges(4,1) + t74 * t46 + t75 * r_base(1) + t77 * (-t9 * pkin(9) + t52) + t70 * (t10 * t44 + t41 * t68) + t71 * t9 - t73 * t68 + t72 * (t10 * t41 - t44 * t68)) * g(1);
U  = t1;
