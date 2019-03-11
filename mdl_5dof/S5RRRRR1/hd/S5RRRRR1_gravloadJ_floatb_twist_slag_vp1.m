% Calculate Gravitation load on the joints for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:37:05
% EndTime: 2019-03-08 18:37:06
% DurationCPUTime: 0.50s
% Computational Cost: add. (360->92), mult. (339->126), div. (0->0), fcn. (292->10), ass. (0->57)
t26 = qJ(2) + qJ(3);
t24 = qJ(4) + t26;
t19 = sin(t24);
t20 = cos(t24);
t27 = sin(qJ(5));
t64 = rSges(6,2) * t27;
t66 = rSges(6,3) + pkin(6);
t77 = t19 * t64 + t20 * t66;
t31 = cos(qJ(2));
t25 = t31 * pkin(2);
t22 = sin(t26);
t23 = cos(t26);
t49 = -t23 * rSges(4,1) + t22 * rSges(4,2);
t76 = -t25 + t49;
t52 = t66 * t19;
t68 = t20 * pkin(4);
t75 = t68 + t52;
t30 = cos(qJ(5));
t65 = rSges(6,1) * t30;
t51 = -pkin(4) - t65;
t44 = t20 * rSges(5,1) - t19 * rSges(5,2);
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t73 = g(1) * t32 + g(2) * t29;
t72 = pkin(3) * t22;
t71 = pkin(3) * t23;
t28 = sin(qJ(2));
t67 = t28 * pkin(2);
t59 = t29 * t27;
t58 = t29 * t30;
t57 = t32 * t27;
t56 = t32 * t30;
t55 = t77 * t29;
t54 = t77 * t32;
t50 = t25 + t71;
t47 = -t31 * rSges(3,1) + t28 * rSges(3,2);
t45 = -rSges(4,1) * t22 - rSges(4,2) * t23;
t43 = -rSges(5,1) * t19 - rSges(5,2) * t20;
t42 = pkin(1) - t47;
t41 = -t44 - t71;
t40 = pkin(1) - t76;
t39 = t43 * t29;
t38 = t43 * t32;
t13 = t20 * t64;
t37 = -t20 * t65 + t13 - t68;
t35 = g(1) * t54 + g(2) * t55;
t34 = (-g(3) * t66 + t73 * t51) * t19;
t11 = -t67 - t72;
t10 = pkin(1) + t50;
t7 = t32 * t11;
t6 = t29 * t11;
t5 = t32 * t10;
t4 = t20 * t56 - t59;
t3 = -t20 * t57 - t58;
t2 = -t20 * t58 - t57;
t1 = t20 * t59 - t56;
t8 = [-m(2) * (g(1) * (-t29 * rSges(2,1) - t32 * rSges(2,2)) + g(2) * (t32 * rSges(2,1) - t29 * rSges(2,2))) - m(3) * ((-g(1) * rSges(3,3) + g(2) * t42) * t32 + (-g(2) * rSges(3,3) - g(1) * t42) * t29) - m(4) * ((-g(1) * rSges(4,3) + g(2) * t40) * t32 + (-g(2) * rSges(4,3) - g(1) * t40) * t29) - m(5) * (g(2) * t5 + (-g(1) * rSges(5,3) + g(2) * t44) * t32 + (g(1) * (-t10 - t44) - g(2) * rSges(5,3)) * t29) - m(6) * ((t4 * rSges(6,1) + t3 * rSges(6,2) + t75 * t32 + t5) * g(2) + (t2 * rSges(6,1) + t1 * rSges(6,2) + (-t10 - t75) * t29) * g(1)) -m(3) * (g(3) * t47 + t73 * (-rSges(3,1) * t28 - rSges(3,2) * t31)) - m(4) * (g(3) * t76 + t73 * (t45 - t67)) - m(5) * (g(1) * (t7 + t38) + g(2) * (t6 + t39) + g(3) * (-t25 + t41)) - m(6) * (g(1) * (t7 + t54) + g(2) * (t6 + t55) + g(3) * (t37 - t50) + t34) -m(6) * t35 + (-m(4) * t49 - m(5) * t41 - m(6) * (t51 * t20 + t13 - t52 - t71)) * g(3) + t73 * (-m(4) * t45 - m(5) * (t43 - t72) - m(6) * (t51 * t19 - t72)) -m(5) * (g(1) * t38 + g(2) * t39 - g(3) * t44) - m(6) * (g(3) * t37 + t34 + t35) -m(6) * (g(1) * (t3 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-t1 * rSges(6,1) + t2 * rSges(6,2)) + g(3) * (rSges(6,1) * t27 + rSges(6,2) * t30) * t19)];
taug  = t8(:);
