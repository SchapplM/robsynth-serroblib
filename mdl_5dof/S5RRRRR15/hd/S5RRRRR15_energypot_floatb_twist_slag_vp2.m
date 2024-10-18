% Calculate potential energy for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_base [3x1]
%   Base position in world frame
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR15_energypot_floatb_twist_slag_vp2(qJ, r_base, g, ...
  pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp2: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:23:52
% EndTime: 2024-09-27 22:23:52
% DurationCPUTime: 0.25s
% Computational Cost: add. (334->122), mult. (372->153), div. (0->0), fcn. (365->26), ass. (0->71)
t36 = qJ(2) + qJ(3);
t29 = pkin(5) + t36;
t23 = qJ(4) + t29;
t14 = sin(t23) / 0.2e1;
t30 = pkin(5) - t36;
t24 = -qJ(4) + t30;
t15 = cos(t24) / 0.2e1;
t51 = pkin(8) + pkin(9);
t35 = pkin(10) + t51;
t39 = cos(pkin(6));
t16 = pkin(11) * t39 + t35;
t17 = sin(t24);
t18 = cos(t23);
t19 = sin(t29);
t20 = sin(t30);
t21 = cos(t29);
t22 = cos(t30);
t38 = sin(pkin(5));
t40 = cos(pkin(5));
t43 = sin(qJ(3));
t44 = sin(qJ(2));
t48 = cos(qJ(3));
t49 = cos(qJ(2));
t52 = pkin(3) * t43 * t49 + (pkin(3) * t48 + pkin(2)) * t44;
t42 = sin(qJ(4));
t47 = cos(qJ(4));
t37 = sin(pkin(6));
t76 = pkin(11) * t37;
t12 = -pkin(4) * t42 + t47 * t76;
t9 = pkin(4) * t47 + t42 * t76 + pkin(3);
t3 = t12 * t43 + t48 * t9 + pkin(2);
t4 = t12 * t48 - t43 * t9;
t53 = t3 * t44 - t4 * t49;
t77 = pkin(2) * t44;
t87 = -(t19 - t20) * mrSges(4,1) / 0.2e1 - (t22 + t21) * mrSges(4,2) / 0.2e1 - m(4) * (-t38 * t51 + t40 * t77) - m(5) * (-t35 * t38 + t40 * t52) + m(6) * (t16 * t38 - t40 * t53) - (t14 - t17 / 0.2e1) * mrSges(5,1) - (t15 + t18 / 0.2e1) * mrSges(5,2) - mrSges(2,2);
t86 = t37 * mrSges(6,3);
t85 = -m(2) - m(3) - m(4) - m(5) - m(6);
t84 = -m(1) + t85;
t82 = m(3) * pkin(8) + t39 * mrSges(6,3) + mrSges(3,3) + mrSges(4,3) + mrSges(5,3);
t28 = t49 * pkin(2) + pkin(1);
t32 = cos(t36);
t80 = -m(3) * pkin(1) - m(4) * t28 - m(5) * (pkin(3) * t32 + t28) - m(6) * (t3 * t49 + t4 * t44 + pkin(1)) - t32 * mrSges(4,1) + sin(t36) * mrSges(4,2) - mrSges(2,1);
t33 = qJ(4) + t36;
t25 = sin(t33);
t45 = sin(qJ(1));
t75 = t25 * t45;
t50 = cos(qJ(1));
t74 = t25 * t50;
t26 = cos(t33);
t73 = t26 * t39;
t72 = t26 * t45;
t71 = t26 * t50;
t70 = t38 * t45;
t69 = t38 * t50;
t41 = sin(qJ(5));
t68 = t41 * t45;
t67 = t41 * t50;
t66 = t44 * t45;
t65 = t44 * t50;
t46 = cos(qJ(5));
t64 = t45 * t46;
t63 = t45 * t49;
t62 = t46 * t50;
t61 = t49 * t50;
t59 = t40 * t68;
t58 = t40 * t67;
t57 = t40 * t64;
t56 = t40 * t62;
t55 = t37 * t70;
t54 = t37 * t69;
t1 = (-m(1) * r_base(3) - mrSges(1,3) - mrSges(2,3) - (t22 / 0.2e1 - t21 / 0.2e1) * mrSges(4,1) - (t19 / 0.2e1 + t20 / 0.2e1) * mrSges(4,2) - (t15 - t18 / 0.2e1) * mrSges(5,1) - (t14 + t17 / 0.2e1) * mrSges(5,2) + (-t44 * mrSges(3,1) - t49 * mrSges(3,2) - m(4) * t77 - m(5) * t52 - m(6) * t53 - (t25 * t46 + t41 * t73) * mrSges(6,1) - (-t25 * t41 + t46 * t73) * mrSges(6,2) + t26 * t86) * t38 + t85 * (pkin(7) + r_base(3)) + (-m(4) * t51 - m(5) * t35 - m(6) * t16 - (mrSges(6,1) * t41 + mrSges(6,2) * t46) * t37 - t82) * t40) * g(3) + (-mrSges(1,2) - (t40 * t65 + t63) * mrSges(3,1) - (t40 * t61 - t66) * mrSges(3,2) - t72 * mrSges(5,1) + t75 * mrSges(5,2) - ((t39 * t58 + t64) * t26 + (-t39 * t68 + t56) * t25 - t41 * t54) * mrSges(6,1) - ((t39 * t56 - t68) * t26 + (-t39 * t64 - t58) * t25 - t46 * t54) * mrSges(6,2) - (-t40 * t71 + t75) * t86 + t80 * t45 + t84 * r_base(2) + t82 * t69 + t87 * t50) * g(2) + (-mrSges(1,1) - (-t40 * t66 + t61) * mrSges(3,1) - (-t40 * t63 - t65) * mrSges(3,2) - t71 * mrSges(5,1) + t74 * mrSges(5,2) - ((-t39 * t59 + t62) * t26 + (-t39 * t67 - t57) * t25 + t41 * t55) * mrSges(6,1) - ((-t39 * t57 - t67) * t26 + (-t39 * t62 + t59) * t25 + t46 * t55) * mrSges(6,2) - (t40 * t72 + t74) * t86 + t80 * t50 + t84 * r_base(1) - t82 * t70 - t87 * t45) * g(1);
U = t1;
