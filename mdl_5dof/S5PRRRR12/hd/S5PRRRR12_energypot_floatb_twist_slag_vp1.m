% Calculate potential energy for
% S5PRRRR12
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
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
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S5PRRRR12_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR12_energypot_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:08
% EndTime: 2024-09-28 18:07:09
% DurationCPUTime: 0.55s
% Computational Cost: add. (334->154), mult. (398->208), div. (0->0), fcn. (399->26), ass. (0->72)
t56 = pkin(7) + pkin(8);
t40 = pkin(9) + t56;
t91 = t40 + rSges(5,3);
t90 = t56 + rSges(4,3);
t42 = sin(pkin(11));
t45 = cos(pkin(11));
t89 = g(1) * t42 - g(2) * t45;
t88 = rSges(3,3) + pkin(7);
t46 = cos(pkin(6));
t87 = t40 + (rSges(6,3) + pkin(10)) * t46;
t43 = sin(pkin(6));
t85 = pkin(10) * t43;
t55 = cos(qJ(2));
t31 = t55 * pkin(2) + pkin(1);
t47 = cos(pkin(5));
t81 = t42 * t47;
t80 = t45 * t47;
t48 = sin(qJ(5));
t78 = t46 * t48;
t52 = cos(qJ(5));
t77 = t46 * t52;
t76 = t47 * t48;
t51 = sin(qJ(2));
t75 = t47 * t51;
t74 = t47 * t52;
t73 = t47 * t55;
t41 = qJ(2) + qJ(3);
t72 = t42 * t85;
t71 = t45 * t85;
t70 = t42 * pkin(1) + r_base(2);
t69 = t45 * pkin(1) + r_base(1);
t68 = t46 * t76;
t67 = t46 * t74;
t66 = qJ(1) + r_base(3);
t33 = pkin(5) - t41;
t32 = pkin(5) + t41;
t64 = rSges(6,1) * t48 + rSges(6,2) * t52;
t50 = sin(qJ(3));
t54 = cos(qJ(3));
t63 = pkin(3) * t50 * t55 + (pkin(3) * t54 + pkin(2)) * t51;
t35 = cos(t41);
t62 = rSges(4,1) * t35 - rSges(4,2) * sin(t41) + t31;
t38 = qJ(4) + t41;
t28 = sin(t38);
t29 = cos(t38);
t61 = rSges(5,1) * t29 - rSges(5,2) * t28 + pkin(3) * t35 + t31;
t49 = sin(qJ(4));
t53 = cos(qJ(4));
t7 = -pkin(4) * t42 + t47 * t71;
t9 = pkin(4) * t80 + t72;
t60 = pkin(3) * t42 + t49 * t9 - t53 * t7;
t6 = pkin(4) * t45 + t47 * t72;
t8 = pkin(4) * t81 - t71;
t59 = pkin(3) * t45 - t49 * t8 + t53 * t6;
t26 = qJ(4) + t32;
t16 = sin(t26) / 0.2e1;
t27 = -qJ(4) + t33;
t17 = cos(t27) / 0.2e1;
t22 = sin(t27);
t23 = cos(t26);
t44 = sin(pkin(5));
t58 = rSges(5,1) * (t16 - t22 / 0.2e1) + rSges(5,2) * (t17 + t23 / 0.2e1) + t47 * t63 - t91 * t44;
t19 = sin(t32) / 0.2e1;
t20 = cos(t33) / 0.2e1;
t24 = sin(t33);
t25 = cos(t32);
t57 = rSges(4,1) * (t19 - t24 / 0.2e1) + rSges(4,2) * (t20 + t25 / 0.2e1) + pkin(2) * t75 - t90 * t44;
t14 = -pkin(4) * t49 + t53 * t85;
t13 = pkin(4) * t53 + t49 * t85 + pkin(3);
t2 = pkin(3) * t80 + t49 * t7 + t53 * t9;
t1 = -pkin(3) * t81 - t49 * t6 - t53 * t8;
t3 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (rSges(2,1) * t45 - rSges(2,2) * t42 + r_base(1)) + g(2) * (rSges(2,1) * t42 + rSges(2,2) * t45 + r_base(2)) + g(3) * (rSges(2,3) + t66)) - m(3) * (g(1) * ((-t42 * t75 + t45 * t55) * rSges(3,1) + (-t42 * t73 - t45 * t51) * rSges(3,2) + t69) + g(2) * ((t42 * t55 + t45 * t75) * rSges(3,1) + (-t42 * t51 + t45 * t73) * rSges(3,2) + t70) + g(3) * (t88 * t47 + t66) + (g(3) * (rSges(3,1) * t51 + rSges(3,2) * t55) + t89 * t88) * t44) - m(4) * (g(1) * (-t42 * t57 + t45 * t62 + r_base(1)) + g(2) * (t42 * t62 + t45 * t57 + r_base(2)) + g(3) * (t44 * t51 * pkin(2) + (t20 - t25 / 0.2e1) * rSges(4,1) + (t19 + t24 / 0.2e1) * rSges(4,2) + t90 * t47 + t66)) - m(5) * (g(1) * (-t42 * t58 + t45 * t61 + r_base(1)) + g(2) * (t42 * t61 + t45 * t58 + r_base(2)) + g(3) * ((t17 - t23 / 0.2e1) * rSges(5,1) + (t16 + t22 / 0.2e1) * rSges(5,2) + t91 * t47 + t63 * t44 + t66)) - m(6) * (g(1) * ((t45 * pkin(2) + t1 * t50 + t54 * t59) * t55 + (-pkin(2) * t81 + t1 * t54 - t50 * t59) * t51 + ((-t42 * t68 + t45 * t52) * t29 + (-t42 * t74 - t45 * t78) * t28) * rSges(6,1) + ((-t42 * t67 - t45 * t48) * t29 + (t42 * t76 - t45 * t77) * t28) * rSges(6,2) + t69) + g(2) * ((t42 * pkin(2) + t2 * t50 + t54 * t60) * t55 + (pkin(2) * t80 + t2 * t54 - t50 * t60) * t51 + ((t42 * t52 + t45 * t68) * t29 + (-t42 * t78 + t45 * t74) * t28) * rSges(6,1) + ((-t42 * t48 + t45 * t67) * t29 + (-t42 * t77 - t45 * t76) * t28) * rSges(6,2) + t70) + g(3) * (t87 * t47 + t66) + (g(3) * t64 * t47 + (g(1) * (t28 * t45 + t29 * t81) + g(2) * (t28 * t42 - t29 * t80)) * rSges(6,3)) * t43 + (g(3) * ((t13 * t54 + t14 * t50 + pkin(2)) * t51 - (-t13 * t50 + t14 * t54) * t55 + (t28 * t52 + t29 * t78) * rSges(6,1) + (-t28 * t48 + t29 * t77) * rSges(6,2) - t29 * t43 * rSges(6,3)) + t89 * (t43 * t64 + t87)) * t44);
U = t3;
