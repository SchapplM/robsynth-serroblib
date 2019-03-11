% Calculate Gravitation load on the joints for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP9_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP9_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:10
% EndTime: 2019-03-09 06:26:12
% DurationCPUTime: 0.75s
% Computational Cost: add. (334->133), mult. (474->177), div. (0->0), fcn. (446->8), ass. (0->56)
t33 = sin(qJ(3));
t36 = cos(qJ(3));
t58 = rSges(5,3) + pkin(8);
t73 = t58 * t36;
t76 = t33 * pkin(3) - t73;
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t48 = g(1) * t34 - g(2) * t37;
t38 = -pkin(9) - pkin(8);
t51 = rSges(7,3) + qJ(6) - t38;
t75 = t51 * t36;
t52 = rSges(6,3) - t38;
t74 = t52 * t36;
t69 = -pkin(1) - pkin(7);
t31 = qJ(4) + qJ(5);
t23 = sin(t31);
t24 = cos(t31);
t56 = t33 * t34;
t10 = t23 * t37 + t24 * t56;
t9 = -t23 * t56 + t24 * t37;
t68 = t9 * rSges(7,1) - t10 * rSges(7,2);
t67 = t9 * rSges(6,1) - t10 * rSges(6,2);
t55 = t33 * t37;
t11 = t23 * t55 + t24 * t34;
t12 = -t23 * t34 + t24 * t55;
t66 = t11 * rSges(7,1) + t12 * rSges(7,2);
t65 = t11 * rSges(6,1) + t12 * rSges(6,2);
t32 = sin(qJ(4));
t64 = pkin(4) * t32;
t63 = pkin(5) * t23;
t60 = g(3) * t36;
t57 = rSges(4,2) * t36;
t35 = cos(qJ(4));
t54 = t34 * t35;
t53 = t35 * t37;
t27 = t35 * pkin(4);
t19 = pkin(5) * t24 + t27;
t50 = t37 * pkin(1) + t34 * qJ(2);
t49 = t37 * pkin(7) + t50;
t46 = rSges(4,1) * t33 + t57;
t45 = -rSges(6,1) * t23 - rSges(6,2) * t24;
t44 = -rSges(7,1) * t23 - rSges(7,2) * t24;
t43 = rSges(5,1) * t35 - rSges(5,2) * t32 + pkin(3);
t15 = t32 * t55 + t54;
t13 = -t32 * t56 + t53;
t22 = t27 + pkin(3);
t42 = rSges(6,1) * t24 - rSges(6,2) * t23 + t22;
t17 = pkin(3) + t19;
t41 = rSges(7,1) * t24 - rSges(7,2) * t23 + t17;
t40 = t33 * t22 - t74;
t39 = t33 * t17 - t75;
t26 = t37 * qJ(2);
t18 = t63 + t64;
t16 = -t32 * t34 + t33 * t53;
t14 = t32 * t37 + t33 * t54;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t34 - rSges(2,2) * t37) + g(2) * (rSges(2,1) * t37 - rSges(2,2) * t34)) - m(3) * (g(1) * (rSges(3,3) * t37 + t26 + (rSges(3,2) - pkin(1)) * t34) + g(2) * (-rSges(3,2) * t37 + rSges(3,3) * t34 + t50)) - m(4) * (g(1) * (rSges(4,1) * t55 + t37 * t57 + t26) + g(2) * (rSges(4,3) * t37 + t49) + (g(1) * (-rSges(4,3) + t69) + g(2) * t46) * t34) - m(5) * ((rSges(5,1) * t14 + rSges(5,2) * t13 + t76 * t34 + t49) * g(2) + (rSges(5,1) * t16 - rSges(5,2) * t15 + t69 * t34 + t76 * t37 + t26) * g(1)) - m(6) * (g(1) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t26) + g(2) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t49) + (g(1) * t40 + g(2) * t64) * t37 + (g(1) * (-t64 + t69) + g(2) * t40) * t34) - m(7) * (g(1) * (rSges(7,1) * t12 - rSges(7,2) * t11 + t26) + g(2) * (rSges(7,1) * t10 + rSges(7,2) * t9 + t49) + (g(1) * t39 + g(2) * t18) * t37 + (g(1) * (-t18 + t69) + g(2) * t39) * t34) (-m(3) - m(4) - m(5) - m(6) - m(7)) * t48, -m(4) * (-g(3) * t46 + t48 * (rSges(4,1) * t36 - rSges(4,2) * t33)) - m(5) * (g(3) * (-t43 * t33 + t73) + t48 * (t58 * t33 + t43 * t36)) - m(6) * (g(3) * (-t42 * t33 + t74) + t48 * (t52 * t33 + t42 * t36)) - m(7) * (g(3) * (-t41 * t33 + t75) + t48 * (t51 * t33 + t41 * t36)) -m(5) * (g(1) * (rSges(5,1) * t13 - rSges(5,2) * t14) + g(2) * (rSges(5,1) * t15 + rSges(5,2) * t16)) - m(6) * (g(1) * (t13 * pkin(4) + t67) + g(2) * (t15 * pkin(4) + t65)) - m(7) * (g(1) * (-t18 * t56 + t19 * t37 + t68) + g(2) * (t18 * t55 + t19 * t34 + t66)) + (-m(5) * (-rSges(5,1) * t32 - rSges(5,2) * t35) - m(6) * (t45 - t64) - m(7) * (-t18 + t44)) * t60, -m(6) * (g(1) * t67 + g(2) * t65) - m(7) * (g(1) * (t9 * pkin(5) + t68) + g(2) * (t11 * pkin(5) + t66)) + (-m(6) * t45 - m(7) * (t44 - t63)) * t60, -m(7) * (g(3) * t33 - t48 * t36)];
taug  = t1(:);
