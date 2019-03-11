% Calculate minimal parameter regressor of gravitation load for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t53 = sin(pkin(13));
t54 = sin(pkin(12));
t44 = t54 * t53;
t57 = cos(pkin(13));
t58 = cos(pkin(12));
t51 = t58 * t57;
t60 = cos(pkin(6));
t37 = -t60 * t51 + t44;
t55 = sin(pkin(7));
t56 = sin(pkin(6));
t48 = t56 * t55;
t59 = cos(pkin(7));
t66 = t37 * t59 + t58 * t48;
t45 = t54 * t57;
t49 = t58 * t53;
t38 = t60 * t45 + t49;
t47 = t56 * t54;
t65 = t38 * t59 - t55 * t47;
t64 = t57 * t59 * t56 + t60 * t55;
t63 = cos(qJ(3));
t29 = qJ(4) + qJ(5);
t28 = cos(t29);
t30 = sin(qJ(6));
t62 = t28 * t30;
t33 = cos(qJ(6));
t61 = t28 * t33;
t50 = t58 * t56;
t46 = t56 * t53;
t22 = t60 * t49 + t45;
t32 = sin(qJ(3));
t12 = t22 * t63 - t66 * t32;
t23 = -t60 * t44 + t51;
t14 = t23 * t63 - t65 * t32;
t16 = t64 * t32 + t63 * t46;
t17 = t37 * t55 - t59 * t50;
t18 = t38 * t55 + t59 * t47;
t21 = -t57 * t48 + t60 * t59;
t27 = sin(t29);
t43 = g(1) * (-t14 * t27 + t18 * t28) + g(2) * (-t12 * t27 + t17 * t28) + g(3) * (-t16 * t27 + t21 * t28);
t11 = t22 * t32 + t66 * t63;
t13 = t23 * t32 + t65 * t63;
t15 = t32 * t46 - t64 * t63;
t42 = g(1) * t13 + g(2) * t11 + g(3) * t15;
t34 = cos(qJ(4));
t31 = sin(qJ(4));
t10 = t16 * t28 + t21 * t27;
t8 = t14 * t28 + t18 * t27;
t6 = t12 * t28 + t17 * t27;
t4 = g(1) * t8 + g(2) * t6 + g(3) * t10;
t2 = t43 * t33;
t1 = t43 * t30;
t3 = [-g(3), -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, -g(1) * t47 + g(2) * t50 - g(3) * t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, t42, g(1) * t14 + g(2) * t12 + g(3) * t16, 0, 0, 0, 0, 0, t42 * t34, -t42 * t31, 0, 0, 0, 0, 0, t42 * t28, -t42 * t27, 0, 0, 0, 0, 0, -g(1) * (-t13 * t61 + t14 * t30) - g(2) * (-t11 * t61 + t12 * t30) - g(3) * (-t15 * t61 + t16 * t30) -g(1) * (t13 * t62 + t14 * t33) - g(2) * (t11 * t62 + t12 * t33) - g(3) * (t15 * t62 + t16 * t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t31 + t18 * t34) - g(2) * (-t12 * t31 + t17 * t34) - g(3) * (-t16 * t31 + t21 * t34) -g(1) * (-t14 * t34 - t18 * t31) - g(2) * (-t12 * t34 - t17 * t31) - g(3) * (-t16 * t34 - t21 * t31) 0, 0, 0, 0, 0, -t43, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t33 - t8 * t30) - g(2) * (t11 * t33 - t6 * t30) - g(3) * (-t10 * t30 + t15 * t33) -g(1) * (-t13 * t30 - t8 * t33) - g(2) * (-t11 * t30 - t6 * t33) - g(3) * (-t10 * t33 - t15 * t30);];
taug_reg  = t3;
