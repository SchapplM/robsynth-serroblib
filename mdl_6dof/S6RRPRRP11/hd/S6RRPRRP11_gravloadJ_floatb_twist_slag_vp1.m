% Calculate Gravitation load on the joints for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRP11_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP11_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:27
% EndTime: 2019-03-09 12:45:28
% DurationCPUTime: 0.75s
% Computational Cost: add. (367->159), mult. (560->217), div. (0->0), fcn. (529->8), ass. (0->66)
t77 = rSges(5,3) + pkin(8);
t45 = -pkin(9) - pkin(8);
t67 = rSges(6,3) - t45;
t66 = rSges(7,3) + qJ(6) - t45;
t44 = cos(qJ(1));
t41 = sin(qJ(1));
t80 = g(2) * t41;
t53 = g(1) * t44 + t80;
t38 = qJ(4) + qJ(5);
t29 = sin(t38);
t30 = cos(t38);
t40 = sin(qJ(2));
t71 = t40 * t44;
t10 = t29 * t71 + t30 * t41;
t9 = -t29 * t41 + t30 * t71;
t84 = t9 * rSges(7,1) - t10 * rSges(7,2);
t72 = t40 * t41;
t11 = t29 * t44 + t30 * t72;
t12 = -t29 * t72 + t30 * t44;
t83 = t11 * rSges(7,1) + t12 * rSges(7,2);
t39 = sin(qJ(4));
t82 = pkin(4) * t39;
t43 = cos(qJ(2));
t79 = g(3) * t43;
t34 = t43 * pkin(2);
t76 = rSges(4,2) * t43;
t19 = pkin(5) * t29 + t82;
t75 = t19 * t40;
t74 = t29 * t43;
t73 = t30 * t43;
t42 = cos(qJ(4));
t70 = t41 * t42;
t69 = t42 * t44;
t68 = t43 * t44;
t33 = t42 * pkin(4);
t20 = pkin(5) * t30 + t33;
t31 = t40 * qJ(3);
t65 = t31 + t34;
t64 = t44 * pkin(1) + t41 * pkin(7);
t63 = qJ(3) * t43;
t62 = -pkin(2) - t77;
t61 = t40 * t82;
t60 = -pkin(2) - t67;
t59 = -pkin(2) - t66;
t58 = -pkin(1) - t31;
t57 = pkin(2) * t68 + t44 * t31 + t64;
t56 = g(1) * t62;
t55 = g(1) * t60;
t54 = g(1) * t59;
t52 = rSges(3,1) * t43 - rSges(3,2) * t40;
t50 = rSges(5,1) * t39 + rSges(5,2) * t42;
t14 = -t39 * t41 + t40 * t69;
t16 = t39 * t44 + t40 * t70;
t49 = rSges(7,1) * t29 + rSges(7,2) * t30 + t19;
t48 = rSges(6,1) * t29 + rSges(6,2) * t30 + t82;
t23 = t41 * t63;
t25 = t44 * t63;
t47 = g(1) * t25 + g(2) * t23 + g(3) * t65;
t46 = g(1) * (t9 * rSges(6,1) - t10 * rSges(6,2)) + g(2) * (t11 * rSges(6,1) + t12 * rSges(6,2)) + g(3) * (-rSges(6,1) * t73 + rSges(6,2) * t74);
t35 = t44 * pkin(7);
t28 = t33 + pkin(3);
t21 = rSges(7,2) * t74;
t18 = pkin(3) + t20;
t17 = -t39 * t72 + t69;
t15 = t39 * t71 + t70;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t41 - rSges(2,2) * t44) + g(2) * (rSges(2,1) * t44 - rSges(2,2) * t41)) - m(3) * (g(1) * (rSges(3,3) * t44 + t35) + g(2) * (rSges(3,1) * t68 - rSges(3,2) * t71 + t64) + (g(1) * (-pkin(1) - t52) + g(2) * rSges(3,3)) * t41) - m(4) * (g(1) * (rSges(4,1) * t44 + t35) + g(2) * (-rSges(4,2) * t68 + rSges(4,3) * t71 + t57) + (g(1) * (-rSges(4,3) * t40 - t34 + t58 + t76) + g(2) * rSges(4,1)) * t41) - m(5) * (g(1) * (rSges(5,1) * t17 - rSges(5,2) * t16 + pkin(3) * t44 + t35) + g(2) * (rSges(5,1) * t15 + rSges(5,2) * t14 + t77 * t68 + t57) + (g(2) * pkin(3) + g(1) * t58 + t43 * t56) * t41) - m(6) * (g(1) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t35) + g(2) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t57) + (g(1) * t28 + g(2) * (t67 * t43 + t61)) * t44 + (g(1) * (t58 - t61) + g(2) * t28 + t43 * t55) * t41) - m(7) * (g(1) * (rSges(7,1) * t12 - rSges(7,2) * t11 + t35) + g(2) * (rSges(7,1) * t10 + rSges(7,2) * t9 + t57) + (g(1) * t18 + g(2) * (t66 * t43 + t75)) * t44 + (g(1) * (t58 - t75) + g(2) * t18 + t43 * t54) * t41) -m(3) * (g(3) * t52 + t53 * (-rSges(3,1) * t40 - rSges(3,2) * t43)) - m(4) * (g(1) * (rSges(4,3) * t68 + t25) + g(2) * (rSges(4,3) * t41 * t43 + t23) + g(3) * (t65 - t76) + (g(3) * rSges(4,3) + t53 * (rSges(4,2) - pkin(2))) * t40) - m(5) * ((g(3) * t77 + t53 * t50) * t43 + (g(3) * t50 + t44 * t56 + t62 * t80) * t40 + t47) - m(6) * ((g(3) * t67 + t53 * t48) * t43 + (g(3) * t48 + t44 * t55 + t60 * t80) * t40 + t47) - m(7) * ((g(3) * t66 + t53 * t49) * t43 + (g(3) * t49 + t44 * t54 + t59 * t80) * t40 + t47) (-m(4) - m(5) - m(6) - m(7)) * (t53 * t40 - t79) -m(5) * (g(1) * (rSges(5,1) * t14 - rSges(5,2) * t15) + g(2) * (rSges(5,1) * t16 + rSges(5,2) * t17) + (-rSges(5,1) * t42 + rSges(5,2) * t39) * t79) - m(6) * ((g(1) * t14 + g(2) * t16 - t42 * t79) * pkin(4) + t46) - m(7) * (g(1) * (-t19 * t41 + t20 * t71 + t84) + g(2) * (t19 * t44 + t20 * t72 + t83) + g(3) * (t21 + (-rSges(7,1) * t30 - t20) * t43)) -m(6) * t46 - m(7) * (g(1) * t84 + g(2) * t83 + g(3) * (-rSges(7,1) * t73 + t21) + (g(1) * t9 + g(2) * t11 - g(3) * t73) * pkin(5)) -m(7) * (g(3) * t40 + t53 * t43)];
taug  = t1(:);
