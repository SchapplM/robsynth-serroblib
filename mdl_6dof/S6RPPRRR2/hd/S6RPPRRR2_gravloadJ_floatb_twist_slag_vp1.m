% Calculate Gravitation load on the joints for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPPRRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:20:19
% EndTime: 2019-03-09 02:20:20
% DurationCPUTime: 0.43s
% Computational Cost: add. (414->101), mult. (328->137), div. (0->0), fcn. (297->12), ass. (0->55)
t66 = rSges(6,3) + pkin(8);
t65 = rSges(7,3) + pkin(9) + pkin(8);
t32 = cos(qJ(5));
t16 = t32 * pkin(5) + pkin(4);
t26 = qJ(5) + qJ(6);
t21 = sin(t26);
t22 = cos(t26);
t30 = sin(qJ(5));
t64 = m(6) * (rSges(6,1) * t32 - rSges(6,2) * t30 + pkin(4)) + m(7) * (rSges(7,1) * t22 - rSges(7,2) * t21 + t16) + m(5) * rSges(5,1);
t24 = pkin(11) + qJ(4);
t19 = cos(t24);
t25 = qJ(1) + pkin(10);
t20 = cos(t25);
t52 = t20 * t22;
t18 = sin(t25);
t57 = t18 * t21;
t5 = t19 * t57 + t52;
t53 = t20 * t21;
t56 = t18 * t22;
t6 = -t19 * t56 + t53;
t63 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t7 = -t19 * t53 + t56;
t8 = t19 * t52 + t57;
t62 = t7 * rSges(7,1) - t8 * rSges(7,2);
t60 = pkin(5) * t30;
t17 = sin(t24);
t59 = g(3) * t17;
t31 = sin(qJ(1));
t58 = t31 * pkin(1);
t55 = t18 * t30;
t54 = t18 * t32;
t51 = t20 * t30;
t50 = t20 * t32;
t29 = -pkin(7) - qJ(3);
t49 = rSges(5,3) - t29;
t28 = cos(pkin(11));
t15 = t28 * pkin(3) + pkin(2);
t33 = cos(qJ(1));
t23 = t33 * pkin(1);
t48 = t20 * t15 + t23;
t47 = rSges(4,3) + qJ(3);
t46 = g(1) * t58;
t45 = -m(4) - m(5) - m(6) - m(7);
t44 = -t29 + t60;
t43 = t19 * rSges(5,1) - t17 * rSges(5,2);
t42 = -rSges(7,1) * t21 - rSges(7,2) * t22;
t41 = rSges(4,1) * t28 - rSges(4,2) * sin(pkin(11)) + pkin(2);
t11 = -t19 * t51 + t54;
t9 = t19 * t55 + t50;
t38 = t19 * pkin(4) + t66 * t17;
t37 = t19 * t16 + t65 * t17;
t36 = m(5) * rSges(5,2) - m(6) * t66 - m(7) * t65;
t12 = t19 * t50 + t55;
t10 = -t19 * t54 + t51;
t1 = [-m(2) * (g(1) * (-t31 * rSges(2,1) - t33 * rSges(2,2)) + g(2) * (t33 * rSges(2,1) - t31 * rSges(2,2))) - m(3) * (g(1) * (-t18 * rSges(3,1) - t20 * rSges(3,2) - t58) + g(2) * (t20 * rSges(3,1) - t18 * rSges(3,2) + t23)) - m(4) * (-t46 + g(2) * t23 + (g(1) * t47 + g(2) * t41) * t20 + (-g(1) * t41 + g(2) * t47) * t18) - m(5) * (-t46 + g(2) * t48 + (g(1) * t49 + g(2) * t43) * t20 + (g(1) * (-t15 - t43) + g(2) * t49) * t18) - m(6) * (g(1) * (t10 * rSges(6,1) + t9 * rSges(6,2) - t58) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t48) + (-g(1) * t29 + g(2) * t38) * t20 + (g(1) * (-t15 - t38) - g(2) * t29) * t18) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2) - t58) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t48) + (g(1) * t44 + g(2) * t37) * t20 + (g(1) * (-t15 - t37) + g(2) * t44) * t18) (-m(3) + t45) * g(3), t45 * (g(1) * t18 - g(2) * t20) (t36 * t17 - t64 * t19) * g(3) + (g(1) * t20 + g(2) * t18) * (t64 * t17 + t36 * t19) -m(6) * (g(1) * (t11 * rSges(6,1) - t12 * rSges(6,2)) + g(2) * (-t9 * rSges(6,1) + t10 * rSges(6,2))) - m(7) * (g(1) * (t11 * pkin(5) + t62) + g(2) * (-t9 * pkin(5) + t63)) + (-m(6) * (-rSges(6,1) * t30 - rSges(6,2) * t32) - m(7) * (t42 - t60)) * t59, -m(7) * (g(1) * t62 + g(2) * t63 + t42 * t59)];
taug  = t1(:);
