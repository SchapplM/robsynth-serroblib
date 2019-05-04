% Calculate Gravitation load on the joints for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(1,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14V3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:27
% EndTime: 2019-04-12 15:03:30
% DurationCPUTime: 0.89s
% Computational Cost: add. (271->156), mult. (694->264), div. (0->0), fcn. (754->10), ass. (0->64)
t31 = sin(qJ(4));
t36 = cos(qJ(4));
t38 = cos(qJ(1));
t53 = t38 * t36;
t33 = sin(qJ(1));
t37 = cos(qJ(2));
t59 = t33 * t37;
t16 = t31 * t59 + t53;
t29 = sin(qJ(6));
t34 = cos(qJ(6));
t54 = t38 * t31;
t57 = t36 * t37;
t17 = t33 * t57 - t54;
t30 = sin(qJ(5));
t35 = cos(qJ(5));
t32 = sin(qJ(2));
t63 = t32 * t33;
t4 = t17 * t35 + t30 * t63;
t77 = -t16 * t34 + t4 * t29;
t76 = -t16 * t29 - t4 * t34;
t28 = t32 * qJ(3);
t45 = t37 * rSges(4,1) + t32 * rSges(4,3);
t75 = -t45 - t28;
t74 = g(1) * t38 + g(2) * t33;
t71 = rSges(7,3) * t30;
t68 = t29 * t35;
t67 = t30 * t31;
t66 = t31 * t32;
t65 = t31 * t35;
t64 = t31 * t37;
t62 = t32 * t35;
t61 = t32 * t36;
t60 = t32 * t38;
t58 = t34 * t35;
t56 = t37 * t30;
t55 = t37 * t35;
t52 = qJ(3) * t37;
t51 = rSges(6,3) * t66;
t50 = t29 * t66;
t49 = t34 * t66;
t48 = t33 * t28;
t47 = t37 * rSges(3,1) - t32 * rSges(3,2);
t44 = -rSges(4,1) * t32 + rSges(4,3) * t37;
t42 = -rSges(6,1) * t35 + rSges(6,2) * t30;
t41 = rSges(7,1) * t34 - rSges(7,2) * t29;
t40 = -t17 * t30 + t33 * t62;
t15 = t35 * t61 - t56;
t39 = t30 * t61 + t55;
t25 = t38 * t52;
t24 = t38 * t28;
t23 = t33 * t52;
t21 = t33 * t31 + t37 * t53;
t20 = -t33 * t36 + t37 * t54;
t19 = t32 * t30 + t36 * t55;
t18 = t36 * t56 - t62;
t12 = t15 * t38;
t11 = t39 * t38;
t10 = t15 * t33;
t9 = t39 * t33;
t8 = t21 * t35 + t30 * t60;
t7 = t21 * t30 - t35 * t60;
t2 = t20 * t29 + t8 * t34;
t1 = t20 * t34 - t8 * t29;
t3 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - t38 * rSges(2,2)) + g(2) * (t38 * rSges(2,1) - t33 * rSges(2,2))) - m(3) * (g(1) * (t38 * rSges(3,3) - t47 * t33) + g(2) * (t33 * rSges(3,3) + t47 * t38)) - m(4) * (g(2) * t24 + (g(1) * rSges(4,2) + g(2) * t45) * t38 + (g(2) * rSges(4,2) + g(1) * t75) * t33) - m(5) * (g(1) * (-t17 * rSges(5,1) + t16 * rSges(5,2) + (-rSges(5,3) - qJ(3)) * t63) + g(2) * (t21 * rSges(5,1) - t20 * rSges(5,2) + rSges(5,3) * t60 + t24)) - m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t40 - t16 * rSges(6,3) - t48) + g(2) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t20 * rSges(6,3) + t24)) - m(7) * (g(1) * (t76 * rSges(7,1) + t77 * rSges(7,2) + t40 * rSges(7,3) - t48) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + t7 * rSges(7,3) + t24)) -m(3) * (g(3) * t47 + t74 * (-rSges(3,1) * t32 - rSges(3,2) * t37)) - m(4) * (g(1) * (t44 * t38 + t25) + g(2) * (t44 * t33 + t23) - g(3) * t75) - m(5) * (g(1) * (t38 * t37 * rSges(5,3) + t25) + g(2) * (rSges(5,3) * t59 + t23) + g(3) * (rSges(5,1) * t57 - rSges(5,2) * t64 + t28) + (g(3) * rSges(5,3) + t74 * (-rSges(5,1) * t36 + rSges(5,2) * t31)) * t32) - m(6) * (g(1) * (-t12 * rSges(6,1) + t11 * rSges(6,2) - t38 * t51 + t25) + g(2) * (-t10 * rSges(6,1) + t9 * rSges(6,2) - t33 * t51 + t23) + g(3) * (t19 * rSges(6,1) - t18 * rSges(6,2) + rSges(6,3) * t64 + t28)) - m(7) * (g(1) * (t25 + (-t12 * t34 - t38 * t50) * rSges(7,1) + (t12 * t29 - t38 * t49) * rSges(7,2) - t11 * rSges(7,3)) + g(2) * (t23 + (-t10 * t34 - t33 * t50) * rSges(7,1) + (t10 * t29 - t33 * t49) * rSges(7,2) - t9 * rSges(7,3)) + g(3) * (t28 + (t19 * t34 + t29 * t64) * rSges(7,1) + (-t19 * t29 + t34 * t64) * rSges(7,2) + t18 * rSges(7,3))) (-m(4) - m(5) - m(6) - m(7)) * (-g(3) * t37 + t74 * t32) -m(5) * (g(1) * (-t20 * rSges(5,1) - t21 * rSges(5,2)) + g(2) * (-t16 * rSges(5,1) - t17 * rSges(5,2))) - m(6) * (g(1) * (t21 * rSges(6,3) + t42 * t20) + g(2) * (t17 * rSges(6,3) + t42 * t16)) - m(7) * (g(1) * ((-t20 * t58 + t21 * t29) * rSges(7,1) + (t20 * t68 + t21 * t34) * rSges(7,2) - t20 * t71) + g(2) * ((-t16 * t58 + t17 * t29) * rSges(7,1) + (t16 * t68 + t17 * t34) * rSges(7,2) - t16 * t71)) + (-m(5) * (-rSges(5,1) * t31 - rSges(5,2) * t36) - m(6) * (-rSges(6,1) * t65 + rSges(6,2) * t67 + rSges(6,3) * t36) - m(7) * ((t29 * t36 - t31 * t58) * rSges(7,1) + (t29 * t65 + t34 * t36) * rSges(7,2) - rSges(7,3) * t67)) * g(3) * t32, -m(6) * (g(1) * (-t7 * rSges(6,1) - t8 * rSges(6,2)) + g(2) * (rSges(6,1) * t40 - t4 * rSges(6,2)) + g(3) * (-rSges(6,1) * t39 - t15 * rSges(6,2))) - m(7) * (g(1) * (t8 * rSges(7,3) - t41 * t7) + g(2) * (t4 * rSges(7,3) + t40 * t41) + g(3) * (t15 * rSges(7,3) - t39 * t41)) -m(7) * (g(1) * (t1 * rSges(7,1) - t2 * rSges(7,2)) + g(2) * (-t77 * rSges(7,1) + t76 * rSges(7,2)) + g(3) * ((-t15 * t29 + t49) * rSges(7,1) + (-t15 * t34 - t50) * rSges(7,2)))];
taug  = t3(:);
