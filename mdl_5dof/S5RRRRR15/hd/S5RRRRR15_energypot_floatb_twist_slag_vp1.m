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
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5RRRRR15_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR15_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:23:52
% EndTime: 2024-09-27 22:23:52
% DurationCPUTime: 0.20s
% Computational Cost: add. (334->120), mult. (360->137), div. (0->0), fcn. (365->26), ass. (0->63)
t51 = pkin(8) + pkin(9);
t35 = pkin(10) + t51;
t37 = sin(pkin(6));
t39 = cos(pkin(6));
t41 = sin(qJ(5));
t46 = cos(qJ(5));
t86 = t35 + (rSges(6,1) * t41 + rSges(6,2) * t46) * t37 + (pkin(11) + rSges(6,3)) * t39;
t36 = qJ(2) + qJ(3);
t33 = qJ(4) + t36;
t26 = cos(t33);
t42 = sin(qJ(4));
t47 = cos(qJ(4));
t78 = pkin(11) * t37;
t12 = -t42 * pkin(4) + t47 * t78;
t43 = sin(qJ(3));
t48 = cos(qJ(3));
t9 = pkin(4) * t47 + t42 * t78 + pkin(3);
t3 = t12 * t43 + t9 * t48 + pkin(2);
t4 = t12 * t48 - t43 * t9;
t44 = sin(qJ(2));
t49 = cos(qJ(2));
t68 = t46 * t26;
t71 = t41 * t26;
t73 = t37 * rSges(6,3);
t25 = sin(t33);
t74 = t25 * t46;
t75 = t25 * t41;
t85 = rSges(6,1) * (t39 * t71 + t74) + rSges(6,2) * (t39 * t68 - t75) - t26 * t73 + t3 * t44 - t4 * t49;
t84 = t35 + rSges(5,3);
t83 = t51 + rSges(4,3);
t38 = sin(pkin(5));
t40 = cos(pkin(5));
t82 = t86 * t38 - t85 * t40;
t81 = rSges(3,3) + pkin(8);
t77 = t44 * pkin(2);
t28 = t49 * pkin(2) + pkin(1);
t45 = sin(qJ(1));
t70 = t45 * t44;
t69 = t45 * t49;
t50 = cos(qJ(1));
t67 = t50 * t44;
t66 = t50 * t49;
t65 = pkin(7) + r_base(3);
t30 = pkin(5) - t36;
t29 = pkin(5) + t36;
t61 = t49 * t43 * pkin(3) + t44 * (t48 * pkin(3) + pkin(2));
t32 = cos(t36);
t58 = t32 * rSges(4,1) - sin(t36) * rSges(4,2) + t28;
t57 = t26 * rSges(5,1) - t25 * rSges(5,2) + pkin(3) * t32 + t28;
t23 = qJ(4) + t29;
t14 = sin(t23) / 0.2e1;
t24 = -qJ(4) + t30;
t15 = cos(t24) / 0.2e1;
t17 = sin(t24);
t18 = cos(t23);
t56 = (t14 - t17 / 0.2e1) * rSges(5,1) + (t15 + t18 / 0.2e1) * rSges(5,2) + t40 * t61 - t84 * t38;
t19 = sin(t29);
t20 = sin(t30);
t21 = cos(t29);
t22 = cos(t30);
t55 = t40 * t77 + (t19 - t20) * rSges(4,1) / 0.2e1 + (t22 + t21) * rSges(4,2) / 0.2e1 - t83 * t38;
t52 = t3 * t49 + t4 * t44 + pkin(1) + (-t39 * t75 + t68) * rSges(6,1) + (-t39 * t74 - t71) * rSges(6,2) + t25 * t73;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t50 * rSges(2,1) - t45 * rSges(2,2) + r_base(1)) + g(2) * (t45 * rSges(2,1) + t50 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t65)) - m(3) * (g(1) * (t50 * pkin(1) + r_base(1) + (-t40 * t70 + t66) * rSges(3,1) + (-t40 * t69 - t67) * rSges(3,2)) + g(2) * (t45 * pkin(1) + r_base(2) + (t40 * t67 + t69) * rSges(3,1) + (t40 * t66 - t70) * rSges(3,2)) + g(3) * (t81 * t40 + t65) + (g(3) * (rSges(3,1) * t44 + rSges(3,2) * t49) + (g(1) * t45 - g(2) * t50) * t81) * t38) - m(4) * (g(1) * (-t45 * t55 + t50 * t58 + r_base(1)) + g(2) * (t45 * t58 + t50 * t55 + r_base(2)) + g(3) * (t38 * t77 + (t22 / 0.2e1 - t21 / 0.2e1) * rSges(4,1) + (t19 / 0.2e1 + t20 / 0.2e1) * rSges(4,2) + t83 * t40 + t65)) - m(5) * (g(1) * (-t45 * t56 + t50 * t57 + r_base(1)) + g(2) * (t45 * t57 + t50 * t56 + r_base(2)) + g(3) * ((t15 - t18 / 0.2e1) * rSges(5,1) + (t14 + t17 / 0.2e1) * rSges(5,2) + t84 * t40 + t61 * t38 + t65)) - m(6) * (g(1) * (t82 * t45 + t52 * t50 + r_base(1)) + g(2) * (t52 * t45 - t82 * t50 + r_base(2)) + g(3) * (t85 * t38 + t86 * t40 + t65));
U = t1;
