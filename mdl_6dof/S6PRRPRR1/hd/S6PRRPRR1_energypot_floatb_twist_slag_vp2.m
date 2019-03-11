% Calculate potential energy for
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRRPRR1_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:50:55
% EndTime: 2019-03-08 21:50:56
% DurationCPUTime: 0.91s
% Computational Cost: add. (378->126), mult. (485->141), div. (0->0), fcn. (534->14), ass. (0->49)
t71 = -m(1) - m(2);
t70 = -m(6) - m(7);
t69 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t41 = -qJ(4) - pkin(8);
t68 = m(4) * pkin(8) - m(5) * t41 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t42 = sin(qJ(6));
t45 = cos(qJ(6));
t67 = -m(7) * pkin(5) - t45 * mrSges(7,1) + t42 * mrSges(7,2) - mrSges(6,1);
t66 = -t42 * mrSges(7,1) - t45 * mrSges(7,2) - mrSges(6,3) - t68;
t43 = sin(qJ(3));
t65 = pkin(3) * t43;
t46 = cos(qJ(3));
t27 = t46 * pkin(3) + pkin(2);
t37 = sin(pkin(11));
t38 = sin(pkin(6));
t64 = t37 * t38;
t39 = cos(pkin(11));
t63 = t38 * t39;
t62 = t38 * t43;
t44 = sin(qJ(2));
t61 = t38 * t44;
t60 = t38 * t46;
t47 = cos(qJ(2));
t59 = t38 * t47;
t40 = cos(pkin(6));
t58 = t40 * t44;
t57 = t40 * t47;
t36 = qJ(3) + pkin(12);
t56 = t37 * pkin(1) + r_base(2);
t55 = t37 * t62;
t54 = qJ(1) + r_base(3);
t53 = t39 * pkin(1) + pkin(7) * t64 + r_base(1);
t52 = t40 * pkin(7) + t54;
t51 = -pkin(7) * t63 + t56;
t29 = cos(t36);
t19 = pkin(4) * t29 + t27;
t28 = sin(t36);
t20 = pkin(4) * t28 + t65;
t35 = -pkin(9) + t41;
t49 = t19 * t61 + t40 * t20 + t35 * t59 + t52;
t30 = qJ(5) + t36;
t26 = cos(t30);
t25 = sin(t30);
t16 = -t37 * t58 + t39 * t47;
t15 = t37 * t57 + t39 * t44;
t14 = t37 * t47 + t39 * t58;
t13 = t37 * t44 - t39 * t57;
t8 = t25 * t40 + t26 * t61;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t54 - mrSges(2,3) - m(6) * t49 - t8 * mrSges(6,1) + mrSges(6,3) * t59 - m(7) * (pkin(5) * t8 + t49) - (-t42 * t59 + t8 * t45) * mrSges(7,1) - (-t8 * t42 - t45 * t59) * mrSges(7,2) + t69 * (t25 * t61 - t40 * t26) + (-m(3) - m(4) - m(5)) * t52 + (-m(5) * t65 - t43 * mrSges(4,1) - t28 * mrSges(5,1) - t46 * mrSges(4,2) - t29 * mrSges(5,2) - mrSges(3,3)) * t40 + (t68 * t47 + (-m(4) * pkin(2) - m(5) * t27 - t46 * mrSges(4,1) - t29 * mrSges(5,1) + t43 * mrSges(4,2) + t28 * mrSges(5,2) - mrSges(3,1)) * t44) * t38) * g(3) + (-m(4) * (pkin(2) * t14 + t51) - m(3) * t51 - m(5) * (t14 * t27 + (-pkin(7) - t65) * t63 + t56) - mrSges(1,2) - t14 * mrSges(3,1) - t37 * mrSges(2,1) - t39 * mrSges(2,2) - (-t14 * t43 - t39 * t60) * mrSges(4,2) - (t14 * t46 - t39 * t62) * mrSges(4,1) - (t14 * t29 - t28 * t63) * mrSges(5,1) + mrSges(3,3) * t63 - (-t14 * t28 - t29 * t63) * mrSges(5,2) + t71 * r_base(2) + t70 * (t14 * t19 - t13 * t35 + (-pkin(7) - t20) * t63 + t56) + t69 * (t14 * t25 + t26 * t63) + t67 * (t14 * t26 - t25 * t63) + t66 * t13) * g(2) + (-m(5) * (pkin(3) * t55 + t16 * t27 + t53) - (t16 * t46 + t55) * mrSges(4,1) - m(4) * (pkin(2) * t16 + t53) - m(3) * t53 - (-t16 * t28 + t29 * t64) * mrSges(5,2) - (t16 * t29 + t28 * t64) * mrSges(5,1) - mrSges(1,1) - mrSges(3,3) * t64 - t16 * mrSges(3,1) + t37 * mrSges(2,2) - t39 * mrSges(2,1) - (-t16 * t43 + t37 * t60) * mrSges(4,2) + t71 * r_base(1) + t70 * (-t15 * t35 + t16 * t19 + t20 * t64 + t53) + t69 * (t16 * t25 - t26 * t64) + t67 * (t16 * t26 + t25 * t64) + t66 * t15) * g(1);
U  = t1;
