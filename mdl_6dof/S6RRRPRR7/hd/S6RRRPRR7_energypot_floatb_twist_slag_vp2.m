% Calculate potential energy for
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPRR7_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPRR7_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:32:30
% EndTime: 2019-03-09 18:32:31
% DurationCPUTime: 0.89s
% Computational Cost: add. (378->126), mult. (485->137), div. (0->0), fcn. (534->14), ass. (0->49)
t71 = -m(2) - m(1);
t70 = -m(6) - m(7);
t69 = -m(7) * pkin(11) + mrSges(6,2) - mrSges(7,3);
t39 = -qJ(4) - pkin(9);
t68 = m(4) * pkin(9) - m(5) * t39 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t40 = sin(qJ(6));
t44 = cos(qJ(6));
t67 = -m(7) * pkin(5) - t44 * mrSges(7,1) + t40 * mrSges(7,2) - mrSges(6,1);
t66 = -t40 * mrSges(7,1) - t44 * mrSges(7,2) - mrSges(6,3) - t68;
t41 = sin(qJ(3));
t65 = pkin(3) * t41;
t45 = cos(qJ(3));
t27 = t45 * pkin(3) + pkin(2);
t37 = sin(pkin(6));
t42 = sin(qJ(2));
t64 = t37 * t42;
t43 = sin(qJ(1));
t63 = t37 * t43;
t46 = cos(qJ(2));
t62 = t37 * t46;
t47 = cos(qJ(1));
t61 = t37 * t47;
t60 = t42 * t43;
t59 = t42 * t47;
t58 = t43 * t46;
t57 = t46 * t47;
t56 = pkin(7) + r_base(3);
t36 = qJ(3) + pkin(12);
t55 = t43 * pkin(1) + r_base(2);
t54 = t41 * t63;
t38 = cos(pkin(6));
t53 = t38 * pkin(8) + t56;
t52 = t47 * pkin(1) + pkin(8) * t63 + r_base(1);
t51 = -pkin(8) * t61 + t55;
t29 = cos(t36);
t19 = pkin(4) * t29 + t27;
t28 = sin(t36);
t20 = pkin(4) * t28 + t65;
t35 = -pkin(10) + t39;
t50 = t19 * t64 + t38 * t20 + t35 * t62 + t53;
t30 = qJ(5) + t36;
t26 = cos(t30);
t25 = sin(t30);
t16 = -t38 * t60 + t57;
t15 = t38 * t58 + t59;
t14 = t38 * t59 + t58;
t13 = -t38 * t57 + t60;
t8 = t25 * t38 + t26 * t64;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t56 - mrSges(2,3) - m(6) * t50 - t8 * mrSges(6,1) + mrSges(6,3) * t62 - m(7) * (pkin(5) * t8 + t50) - (-t40 * t62 + t44 * t8) * mrSges(7,1) - (-t40 * t8 - t44 * t62) * mrSges(7,2) + t69 * (t25 * t64 - t38 * t26) + (-m(3) - m(4) - m(5)) * t53 + (-m(5) * t65 - t41 * mrSges(4,1) - t28 * mrSges(5,1) - t45 * mrSges(4,2) - t29 * mrSges(5,2) - mrSges(3,3)) * t38 + (t68 * t46 + (-m(4) * pkin(2) - m(5) * t27 - t45 * mrSges(4,1) - t29 * mrSges(5,1) + t41 * mrSges(4,2) + t28 * mrSges(5,2) - mrSges(3,1)) * t42) * t37) * g(3) + (-m(4) * (t14 * pkin(2) + t51) - m(3) * t51 - t43 * mrSges(2,1) - mrSges(1,2) - t47 * mrSges(2,2) + mrSges(3,3) * t61 - (-t14 * t41 - t45 * t61) * mrSges(4,2) - (t14 * t45 - t41 * t61) * mrSges(4,1) - (-t14 * t28 - t29 * t61) * mrSges(5,2) - (t14 * t29 - t28 * t61) * mrSges(5,1) - m(5) * (t14 * t27 + (-pkin(8) - t65) * t61 + t55) - t14 * mrSges(3,1) + t71 * r_base(2) + t70 * (t14 * t19 - t13 * t35 + (-pkin(8) - t20) * t61 + t55) + t69 * (t14 * t25 + t26 * t61) + t67 * (t14 * t26 - t25 * t61) + t66 * t13) * g(2) + (-t16 * mrSges(3,1) + t43 * mrSges(2,2) - mrSges(1,1) - m(4) * (pkin(2) * t16 + t52) - m(3) * t52 - t47 * mrSges(2,1) - m(5) * (pkin(3) * t54 + t16 * t27 + t52) - (t16 * t45 + t54) * mrSges(4,1) - (-t16 * t41 + t45 * t63) * mrSges(4,2) - (-t16 * t28 + t29 * t63) * mrSges(5,2) - (t16 * t29 + t28 * t63) * mrSges(5,1) - mrSges(3,3) * t63 + t71 * r_base(1) + t70 * (-t15 * t35 + t16 * t19 + t20 * t63 + t52) + t69 * (t16 * t25 - t26 * t63) + t67 * (t16 * t26 + t25 * t63) + t66 * t15) * g(1);
U  = t1;
