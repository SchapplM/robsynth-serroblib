% Calculate Gravitation load on the joints for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:44
% EndTime: 2019-03-09 04:07:46
% DurationCPUTime: 0.65s
% Computational Cost: add. (324->113), mult. (409->155), div. (0->0), fcn. (379->10), ass. (0->56)
t33 = sin(qJ(3));
t35 = cos(qJ(3));
t30 = sin(pkin(10));
t31 = cos(pkin(10));
t42 = rSges(5,1) * t31 - rSges(5,2) * t30 + pkin(3);
t51 = rSges(5,3) + qJ(4);
t75 = -t42 * t33 + t51 * t35;
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t47 = g(1) * t34 - g(2) * t36;
t32 = -pkin(8) - qJ(4);
t53 = rSges(7,3) + pkin(9) - t32;
t74 = t53 * t35;
t54 = rSges(6,3) - t32;
t73 = t54 * t35;
t68 = -pkin(1) - pkin(7);
t29 = pkin(10) + qJ(5);
t22 = qJ(6) + t29;
t18 = cos(t22);
t17 = sin(t22);
t58 = t34 * t17;
t5 = t18 * t36 - t33 * t58;
t57 = t34 * t18;
t6 = t17 * t36 + t33 * t57;
t67 = t5 * rSges(7,1) - t6 * rSges(7,2);
t59 = t33 * t36;
t7 = t17 * t59 + t57;
t8 = t18 * t59 - t58;
t66 = t7 * rSges(7,1) + t8 * rSges(7,2);
t65 = pkin(4) * t30;
t20 = sin(t29);
t64 = pkin(5) * t20;
t61 = g(3) * t35;
t19 = t31 * pkin(4) + pkin(3);
t60 = rSges(4,2) * t35;
t56 = t34 * t20;
t21 = cos(t29);
t55 = t34 * t21;
t52 = t36 * pkin(1) + t34 * qJ(2);
t50 = -m(5) - m(6) - m(7);
t49 = t36 * pkin(7) + t52;
t45 = rSges(4,1) * t33 + t60;
t44 = rSges(5,1) * t30 + rSges(5,2) * t31;
t43 = -rSges(7,1) * t17 - rSges(7,2) * t18;
t11 = t20 * t59 + t55;
t9 = t21 * t36 - t33 * t56;
t41 = rSges(6,1) * t21 - rSges(6,2) * t20 + t19;
t14 = pkin(5) * t21 + t19;
t40 = rSges(7,1) * t18 - rSges(7,2) * t17 + t14;
t38 = t33 * t19 - t73;
t37 = t33 * t14 - t74;
t25 = t36 * qJ(2);
t15 = t64 + t65;
t12 = t21 * t59 - t56;
t10 = t20 * t36 + t33 * t55;
t1 = [-m(2) * (g(1) * (-t34 * rSges(2,1) - rSges(2,2) * t36) + g(2) * (rSges(2,1) * t36 - t34 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t36 + t25 + (rSges(3,2) - pkin(1)) * t34) + g(2) * (-rSges(3,2) * t36 + t34 * rSges(3,3) + t52)) - m(4) * (g(1) * (rSges(4,1) * t59 + t36 * t60 + t25) + g(2) * (rSges(4,3) * t36 + t49) + (g(1) * (-rSges(4,3) + t68) + g(2) * t45) * t34) - m(5) * (g(1) * t25 + g(2) * t49 + (-g(1) * t75 + g(2) * t44) * t36 + (g(1) * (-t44 + t68) - t75 * g(2)) * t34) - m(6) * (g(1) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t25) + g(2) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t49) + (g(1) * t38 + g(2) * t65) * t36 + (g(1) * (-t65 + t68) + g(2) * t38) * t34) - m(7) * (g(1) * (t8 * rSges(7,1) - t7 * rSges(7,2) + t25) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t49) + (g(1) * t37 + g(2) * t15) * t36 + (g(1) * (-t15 + t68) + g(2) * t37) * t34) (-m(3) - m(4) + t50) * t47, -m(4) * (-g(3) * t45 + t47 * (rSges(4,1) * t35 - rSges(4,2) * t33)) - m(5) * (g(3) * t75 + t47 * (t51 * t33 + t42 * t35)) - m(6) * (g(3) * (-t41 * t33 + t73) + t47 * (t54 * t33 + t41 * t35)) - m(7) * (g(3) * (-t40 * t33 + t74) + t47 * (t53 * t33 + t40 * t35)) t50 * (g(3) * t33 - t47 * t35) -m(6) * (g(1) * (rSges(6,1) * t9 - rSges(6,2) * t10) + g(2) * (rSges(6,1) * t11 + rSges(6,2) * t12)) - m(7) * (g(1) * (t9 * pkin(5) + t67) + g(2) * (t11 * pkin(5) + t66)) + (-m(6) * (-rSges(6,1) * t20 - rSges(6,2) * t21) - m(7) * (t43 - t64)) * t61, -m(7) * (g(1) * t67 + g(2) * t66 + t43 * t61)];
taug  = t1(:);
