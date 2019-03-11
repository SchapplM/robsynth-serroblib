% Calculate potential energy for
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U = S6PRPRRP2_energypot_floatb_twist_slag_vp1(qJ, r_base, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(r_base) && all(size(r_base) == [3 1]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp1: r_base has to be [3x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP2_energypot_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From energy_potential_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:59:50
% EndTime: 2019-03-08 19:59:50
% DurationCPUTime: 0.51s
% Computational Cost: add. (439->124), mult. (938->155), div. (0->0), fcn. (1177->12), ass. (0->57)
t43 = sin(pkin(11));
t46 = cos(pkin(11));
t51 = sin(qJ(2));
t54 = cos(qJ(2));
t33 = -t51 * t43 + t54 * t46;
t79 = rSges(7,1) + pkin(5);
t78 = rSges(7,2) + pkin(9);
t76 = rSges(5,3) + pkin(8);
t75 = rSges(6,3) + pkin(9);
t44 = sin(pkin(10));
t45 = sin(pkin(6));
t74 = t44 * t45;
t47 = cos(pkin(10));
t73 = t47 * t45;
t48 = cos(pkin(6));
t72 = t48 * t51;
t71 = t48 * t54;
t68 = rSges(7,3) + qJ(6);
t40 = t54 * pkin(2) + pkin(1);
t67 = t47 * t40 + r_base(1);
t66 = qJ(1) + r_base(3);
t31 = pkin(2) * t72 + (-pkin(7) - qJ(3)) * t45;
t65 = t47 * t31 + t44 * t40 + r_base(2);
t64 = t48 * pkin(7) + t66;
t62 = t54 * t43 + t51 * t46;
t30 = t62 * t48;
t18 = t47 * t30 + t44 * t33;
t63 = t18 * pkin(3) + t65;
t61 = t45 * t51 * pkin(2) + t48 * qJ(3) + t64;
t20 = -t44 * t30 + t47 * t33;
t60 = t20 * pkin(3) - t44 * t31 + t67;
t59 = t33 * t48;
t29 = t62 * t45;
t58 = t29 * pkin(3) + t61;
t50 = sin(qJ(4));
t53 = cos(qJ(4));
t10 = t18 * t53 - t50 * t73;
t17 = -t44 * t62 + t47 * t59;
t57 = t10 * pkin(4) - t17 * pkin(8) + t63;
t12 = t20 * t53 + t50 * t74;
t19 = -t44 * t59 - t47 * t62;
t56 = t12 * pkin(4) - t19 * pkin(8) + t60;
t23 = t29 * t53 + t48 * t50;
t28 = t33 * t45;
t55 = t23 * pkin(4) - t28 * pkin(8) + t58;
t52 = cos(qJ(5));
t49 = sin(qJ(5));
t22 = t29 * t50 - t48 * t53;
t11 = t20 * t50 - t53 * t74;
t9 = t18 * t50 + t53 * t73;
t6 = t23 * t52 - t28 * t49;
t5 = t23 * t49 + t28 * t52;
t4 = t12 * t52 - t19 * t49;
t3 = t12 * t49 + t19 * t52;
t2 = t10 * t52 - t17 * t49;
t1 = t10 * t49 + t17 * t52;
t7 = -m(1) * (g(1) * (r_base(1) + rSges(1,1)) + g(2) * (r_base(2) + rSges(1,2)) + g(3) * (r_base(3) + rSges(1,3))) - m(2) * (g(1) * (t47 * rSges(2,1) - t44 * rSges(2,2) + r_base(1)) + g(2) * (t44 * rSges(2,1) + t47 * rSges(2,2) + r_base(2)) + g(3) * (rSges(2,3) + t66)) - m(3) * (g(1) * (t47 * pkin(1) + r_base(1) + (-t44 * t72 + t47 * t54) * rSges(3,1) + (-t44 * t71 - t47 * t51) * rSges(3,2)) + g(2) * (t44 * pkin(1) + r_base(2) + (t44 * t54 + t47 * t72) * rSges(3,1) + (-t44 * t51 + t47 * t71) * rSges(3,2)) + g(3) * (t48 * rSges(3,3) + t64) + (g(3) * (rSges(3,1) * t51 + rSges(3,2) * t54) + (g(1) * t44 - g(2) * t47) * (rSges(3,3) + pkin(7))) * t45) - m(4) * (g(1) * (t20 * rSges(4,1) + t19 * rSges(4,2) + (rSges(4,3) * t45 - t31) * t44 + t67) + g(2) * (t18 * rSges(4,1) + t17 * rSges(4,2) - rSges(4,3) * t73 + t65) + g(3) * (t29 * rSges(4,1) + t28 * rSges(4,2) + t48 * rSges(4,3) + t61)) - m(5) * (g(1) * (t12 * rSges(5,1) - t11 * rSges(5,2) - t76 * t19 + t60) + g(2) * (t10 * rSges(5,1) - t9 * rSges(5,2) - t17 * t76 + t63) + g(3) * (t23 * rSges(5,1) - t22 * rSges(5,2) - t28 * t76 + t58)) - m(6) * (g(1) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t11 * t75 + t56) + g(2) * (t2 * rSges(6,1) - t1 * rSges(6,2) + t75 * t9 + t57) + g(3) * (t6 * rSges(6,1) - t5 * rSges(6,2) + t22 * t75 + t55)) - m(7) * (g(1) * (t78 * t11 + t68 * t3 + t4 * t79 + t56) + g(2) * (t68 * t1 + t2 * t79 + t78 * t9 + t57) + g(3) * (t78 * t22 + t68 * t5 + t6 * t79 + t55));
U  = t7;
