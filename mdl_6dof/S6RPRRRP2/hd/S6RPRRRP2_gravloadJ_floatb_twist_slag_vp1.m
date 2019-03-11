% Calculate Gravitation load on the joints for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:59:03
% EndTime: 2019-03-09 05:59:04
% DurationCPUTime: 0.70s
% Computational Cost: add. (487->133), mult. (462->181), div. (0->0), fcn. (434->10), ass. (0->60)
t38 = cos(qJ(3));
t35 = sin(qJ(3));
t64 = rSges(5,3) + pkin(8);
t54 = t64 * t35;
t79 = t38 * pkin(3) + t54;
t40 = -pkin(9) - pkin(8);
t56 = rSges(7,3) + qJ(6) - t40;
t78 = t56 * t35;
t57 = rSges(6,3) - t40;
t77 = t57 * t35;
t32 = qJ(1) + pkin(10);
t25 = sin(t32);
t26 = cos(t32);
t76 = g(1) * t26 + g(2) * t25;
t33 = qJ(4) + qJ(5);
t27 = sin(t33);
t28 = cos(t33);
t60 = t28 * t38;
t10 = -t25 * t60 + t26 * t27;
t61 = t27 * t38;
t9 = t25 * t61 + t26 * t28;
t75 = -t9 * rSges(6,1) + t10 * rSges(6,2);
t74 = -t9 * rSges(7,1) + t10 * rSges(7,2);
t11 = t25 * t28 - t26 * t61;
t12 = t25 * t27 + t26 * t60;
t73 = t11 * rSges(7,1) - t12 * rSges(7,2);
t72 = t11 * rSges(6,1) - t12 * rSges(6,2);
t71 = pkin(5) * t27;
t68 = g(3) * t35;
t34 = sin(qJ(4));
t67 = t34 * pkin(4);
t36 = sin(qJ(1));
t66 = t36 * pkin(1);
t63 = rSges(4,2) * t35;
t62 = t26 * t38;
t59 = t34 * t38;
t37 = cos(qJ(4));
t58 = t37 * t38;
t29 = t37 * pkin(4);
t19 = pkin(5) * t28 + t29;
t39 = cos(qJ(1));
t30 = t39 * pkin(1);
t55 = t26 * pkin(2) + t25 * pkin(7) + t30;
t53 = t26 * pkin(7) - t66;
t52 = rSges(4,1) * t38 - t63;
t50 = -rSges(6,1) * t27 - rSges(6,2) * t28;
t49 = -rSges(7,1) * t27 - rSges(7,2) * t28;
t48 = rSges(5,1) * t37 - rSges(5,2) * t34 + pkin(3);
t15 = t25 * t37 - t26 * t59;
t13 = t25 * t59 + t26 * t37;
t24 = t29 + pkin(3);
t47 = rSges(6,1) * t28 - rSges(6,2) * t27 + t24;
t17 = pkin(3) + t19;
t46 = rSges(7,1) * t28 - rSges(7,2) * t27 + t17;
t45 = t38 * t24 + t77;
t44 = t17 * t38 + t78;
t18 = t67 + t71;
t16 = t25 * t34 + t26 * t58;
t14 = -t25 * t58 + t26 * t34;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t36 - rSges(2,2) * t39) + g(2) * (rSges(2,1) * t39 - rSges(2,2) * t36)) - m(3) * (g(1) * (-rSges(3,1) * t25 - rSges(3,2) * t26 - t66) + g(2) * (rSges(3,1) * t26 - rSges(3,2) * t25 + t30)) - m(4) * (g(1) * (rSges(4,3) * t26 + t53) + g(2) * (rSges(4,1) * t62 - t26 * t63 + t55) + (g(1) * (-pkin(2) - t52) + g(2) * rSges(4,3)) * t25) - m(5) * ((rSges(5,1) * t16 + rSges(5,2) * t15 + t79 * t26 + t55) * g(2) + (rSges(5,1) * t14 + rSges(5,2) * t13 + t53 + (-pkin(2) - t79) * t25) * g(1)) - m(6) * (g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t53) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t55) + (g(1) * t67 + g(2) * t45) * t26 + (g(1) * (-pkin(2) - t45) + g(2) * t67) * t25) - m(7) * (g(1) * (rSges(7,1) * t10 + rSges(7,2) * t9 + t53) + g(2) * (rSges(7,1) * t12 + rSges(7,2) * t11 + t55) + (g(1) * t18 + g(2) * t44) * t26 + (g(1) * (-pkin(2) - t44) + g(2) * t18) * t25) (-m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(4) * (g(3) * t52 + t76 * (-rSges(4,1) * t35 - rSges(4,2) * t38)) - m(5) * (g(3) * (t48 * t38 + t54) + t76 * (-t48 * t35 + t64 * t38)) - m(6) * (g(3) * (t47 * t38 + t77) + t76 * (-t47 * t35 + t57 * t38)) - m(7) * (g(3) * (t46 * t38 + t78) + t76 * (-t46 * t35 + t56 * t38)) -m(5) * (g(1) * (rSges(5,1) * t15 - rSges(5,2) * t16) + g(2) * (-rSges(5,1) * t13 + rSges(5,2) * t14)) - m(6) * (g(1) * (t15 * pkin(4) + t72) + g(2) * (-t13 * pkin(4) + t75)) - m(7) * (g(1) * (-t18 * t62 + t19 * t25 + t73) + g(2) * (-t18 * t25 * t38 - t19 * t26 + t74)) + (-m(5) * (-rSges(5,1) * t34 - rSges(5,2) * t37) - m(6) * (t50 - t67) - m(7) * (-t18 + t49)) * t68, -m(6) * (g(1) * t72 + g(2) * t75) - m(7) * (g(1) * (t11 * pkin(5) + t73) + g(2) * (-t9 * pkin(5) + t74)) + (-m(6) * t50 - m(7) * (t49 - t71)) * t68, -m(7) * (-g(3) * t38 + t35 * t76)];
taug  = t1(:);
