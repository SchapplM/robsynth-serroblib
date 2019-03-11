% Calculate potential energy for
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRR3_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR3_energypot_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:34
% EndTime: 2019-03-08 20:31:35
% DurationCPUTime: 0.89s
% Computational Cost: add. (378->126), mult. (485->139), div. (0->0), fcn. (534->14), ass. (0->47)
t69 = -m(1) - m(2);
t68 = -m(6) - m(7);
t67 = -m(7) * pkin(10) + mrSges(6,2) - mrSges(7,3);
t43 = -pkin(8) - qJ(3);
t66 = m(4) * qJ(3) - m(5) * t43 - mrSges(3,2) + mrSges(4,3) + mrSges(5,3);
t44 = sin(qJ(6));
t46 = cos(qJ(6));
t65 = -m(7) * pkin(5) - t46 * mrSges(7,1) + t44 * mrSges(7,2) - mrSges(6,1);
t64 = -t44 * mrSges(7,1) - t46 * mrSges(7,2) - mrSges(6,3) - t66;
t37 = sin(pkin(12));
t63 = pkin(3) * t37;
t40 = cos(pkin(12));
t27 = t40 * pkin(3) + pkin(2);
t38 = sin(pkin(11));
t39 = sin(pkin(6));
t62 = t38 * t39;
t41 = cos(pkin(11));
t61 = t39 * t41;
t45 = sin(qJ(2));
t60 = t39 * t45;
t47 = cos(qJ(2));
t59 = t39 * t47;
t42 = cos(pkin(6));
t58 = t42 * t45;
t57 = t42 * t47;
t36 = pkin(12) + qJ(4);
t56 = t38 * pkin(1) + r_base(2);
t55 = t37 * t62;
t54 = qJ(1) + r_base(3);
t53 = t41 * pkin(1) + pkin(7) * t62 + r_base(1);
t52 = t42 * pkin(7) + t54;
t51 = -pkin(7) * t61 + t56;
t29 = cos(t36);
t19 = pkin(4) * t29 + t27;
t28 = sin(t36);
t20 = pkin(4) * t28 + t63;
t35 = -pkin(9) + t43;
t49 = t19 * t60 + t42 * t20 + t35 * t59 + t52;
t30 = qJ(5) + t36;
t26 = cos(t30);
t25 = sin(t30);
t16 = -t38 * t58 + t41 * t47;
t15 = t38 * t57 + t41 * t45;
t14 = t38 * t47 + t41 * t58;
t13 = t38 * t45 - t41 * t57;
t8 = t25 * t42 + t26 * t60;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - m(2) * t54 - mrSges(2,3) - m(6) * t49 - t8 * mrSges(6,1) + mrSges(6,3) * t59 - m(7) * (pkin(5) * t8 + t49) - (-t44 * t59 + t8 * t46) * mrSges(7,1) - (-t8 * t44 - t46 * t59) * mrSges(7,2) + t67 * (t25 * t60 - t42 * t26) + (-m(3) - m(4) - m(5)) * t52 + (-m(5) * t63 - t37 * mrSges(4,1) - t28 * mrSges(5,1) - t40 * mrSges(4,2) - t29 * mrSges(5,2) - mrSges(3,3)) * t42 + (t66 * t47 + (-m(4) * pkin(2) - m(5) * t27 - t40 * mrSges(4,1) - t29 * mrSges(5,1) + t37 * mrSges(4,2) + t28 * mrSges(5,2) - mrSges(3,1)) * t45) * t39) * g(3) + (-m(4) * (pkin(2) * t14 + t51) - m(3) * t51 - mrSges(1,2) + mrSges(3,3) * t61 - (-t14 * t37 - t40 * t61) * mrSges(4,2) - (t14 * t40 - t37 * t61) * mrSges(4,1) - (-t14 * t28 - t29 * t61) * mrSges(5,2) - (t14 * t29 - t28 * t61) * mrSges(5,1) - m(5) * (t14 * t27 + (-pkin(7) - t63) * t61 + t56) - t14 * mrSges(3,1) - t38 * mrSges(2,1) - t41 * mrSges(2,2) + t69 * r_base(2) + t68 * (t14 * t19 - t13 * t35 + (-pkin(7) - t20) * t61 + t56) + t67 * (t14 * t25 + t26 * t61) + t65 * (t14 * t26 - t25 * t61) + t64 * t13) * g(2) + (-m(4) * (pkin(2) * t16 + t53) - m(3) * t53 - mrSges(1,1) - m(5) * (pkin(3) * t55 + t16 * t27 + t53) - (t16 * t40 + t55) * mrSges(4,1) - (-t16 * t37 + t40 * t62) * mrSges(4,2) - (-t16 * t28 + t29 * t62) * mrSges(5,2) - (t16 * t29 + t28 * t62) * mrSges(5,1) - mrSges(3,3) * t62 - t16 * mrSges(3,1) + t38 * mrSges(2,2) - t41 * mrSges(2,1) + t69 * r_base(1) + t68 * (-t15 * t35 + t16 * t19 + t20 * t62 + t53) + t67 * (t16 * t25 - t26 * t62) + t65 * (t16 * t26 + t25 * t62) + t64 * t15) * g(1);
U  = t1;
