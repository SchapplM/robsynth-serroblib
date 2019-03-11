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
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% U [1x1]
%   Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6RRRPPR6_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR6_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:39
% EndTime: 2019-03-09 15:47:40
% DurationCPUTime: 0.41s
% Computational Cost: add. (353->130), mult. (523->157), div. (0->0), fcn. (603->12), ass. (0->49)
t38 = sin(qJ(3));
t67 = pkin(3) * t38;
t66 = pkin(9) + rSges(4,3);
t65 = pkin(10) + rSges(7,3);
t34 = sin(pkin(6));
t39 = sin(qJ(2));
t64 = t34 * t39;
t40 = sin(qJ(1));
t63 = t34 * t40;
t43 = cos(qJ(2));
t62 = t34 * t43;
t44 = cos(qJ(1));
t61 = t34 * t44;
t60 = t40 * t39;
t59 = t40 * t43;
t58 = t44 * t39;
t57 = t44 * t43;
t56 = rSges(6,3) + qJ(5);
t55 = pkin(7) + r_base(3);
t54 = t40 * pkin(1) + r_base(2);
t53 = t38 * t63;
t35 = cos(pkin(6));
t52 = t35 * pkin(8) + t55;
t51 = t44 * pkin(1) + pkin(8) * t63 + r_base(1);
t42 = cos(qJ(3));
t27 = t42 * pkin(3) + pkin(2);
t36 = -qJ(4) - pkin(9);
t50 = t27 * t64 + t35 * t67 + t36 * t62 + t52;
t16 = t35 * t59 + t58;
t17 = -t35 * t60 + t57;
t49 = pkin(3) * t53 - t16 * t36 + t17 * t27 + t51;
t33 = qJ(3) + pkin(11);
t28 = sin(t33);
t29 = cos(t33);
t11 = t35 * t28 + t29 * t64;
t48 = t11 * pkin(4) + t50;
t6 = t17 * t29 + t28 * t63;
t47 = t6 * pkin(4) + t49;
t14 = -t35 * t57 + t60;
t15 = t35 * t58 + t59;
t46 = t15 * t27 + (-pkin(8) - t67) * t61 - t14 * t36 + t54;
t4 = t15 * t29 - t28 * t61;
t45 = t4 * pkin(4) + t46;
t41 = cos(qJ(6));
t37 = sin(qJ(6));
t10 = t28 * t64 - t35 * t29;
t5 = t17 * t28 - t29 * t63;
t3 = t15 * t28 + t29 * t61;
t1 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t44 * rSges(2,1) - t40 * rSges(2,2) + r_base(1)) + g(2) * (t40 * rSges(2,1) + t44 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t55)) - m(3) * (g(1) * (t17 * rSges(3,1) - t16 * rSges(3,2) + t51) + g(2) * (t15 * rSges(3,1) - t14 * rSges(3,2) + t54) + g(3) * (t35 * rSges(3,3) + t52) + (g(1) * rSges(3,3) * t40 + g(3) * (rSges(3,1) * t39 + rSges(3,2) * t43) + g(2) * (-rSges(3,3) - pkin(8)) * t44) * t34) - m(4) * (g(1) * (t17 * pkin(2) + (t17 * t42 + t53) * rSges(4,1) + (-t17 * t38 + t42 * t63) * rSges(4,2) + t66 * t16 + t51) + g(2) * (t15 * pkin(2) - pkin(8) * t61 + (t15 * t42 - t38 * t61) * rSges(4,1) + (-t15 * t38 - t42 * t61) * rSges(4,2) + t66 * t14 + t54) + g(3) * ((t38 * rSges(4,1) + t42 * rSges(4,2)) * t35 + (-t66 * t43 + (t42 * rSges(4,1) - t38 * rSges(4,2) + pkin(2)) * t39) * t34 + t52)) - m(5) * (g(1) * (t6 * rSges(5,1) - t5 * rSges(5,2) + t16 * rSges(5,3) + t49) + g(2) * (t4 * rSges(5,1) - t3 * rSges(5,2) + t14 * rSges(5,3) + t46) + g(3) * (t11 * rSges(5,1) - t10 * rSges(5,2) - rSges(5,3) * t62 + t50)) - m(6) * (g(1) * (t16 * rSges(6,1) - t6 * rSges(6,2) + t56 * t5 + t47) + g(2) * (t14 * rSges(6,1) - t4 * rSges(6,2) + t56 * t3 + t45) + g(3) * (-rSges(6,1) * t62 - t11 * rSges(6,2) + t56 * t10 + t48)) - m(7) * (g(1) * (t16 * pkin(5) + t5 * qJ(5) + (t16 * t41 + t5 * t37) * rSges(7,1) + (-t16 * t37 + t5 * t41) * rSges(7,2) + t65 * t6 + t47) + g(2) * (t14 * pkin(5) + t3 * qJ(5) + (t14 * t41 + t3 * t37) * rSges(7,1) + (-t14 * t37 + t3 * t41) * rSges(7,2) + t65 * t4 + t45) + g(3) * (-pkin(5) * t62 + t10 * qJ(5) + (t10 * t37 - t41 * t62) * rSges(7,1) + (t10 * t41 + t37 * t62) * rSges(7,2) + t65 * t11 + t48));
U  = t1;
