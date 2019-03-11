% Calculate minimal parameter regressor of gravitation load for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPPRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t27 = sin(pkin(11));
t34 = sin(qJ(2));
t37 = cos(qJ(2));
t47 = cos(pkin(11));
t21 = -t37 * t27 - t34 * t47;
t42 = -t34 * t27 + t37 * t47;
t28 = sin(pkin(10));
t30 = cos(pkin(10));
t31 = cos(pkin(6));
t48 = t21 * t31;
t12 = t28 * t48 + t30 * t42;
t7 = -t28 * t42 + t30 * t48;
t59 = t28 * t34;
t29 = sin(pkin(6));
t33 = sin(qJ(5));
t58 = t29 * t33;
t36 = cos(qJ(5));
t57 = t29 * t36;
t56 = t29 * t37;
t54 = t31 * t34;
t53 = t31 * t37;
t32 = sin(qJ(6));
t52 = t32 * t33;
t35 = cos(qJ(6));
t51 = t33 * t35;
t46 = t30 * t53;
t43 = -t28 * t53 - t30 * t34;
t39 = t42 * t31;
t11 = t30 * t21 - t28 * t39;
t18 = t42 * t29;
t8 = t28 * t21 + t30 * t39;
t41 = g(1) * (-t11 * t36 - t28 * t58) + g(2) * (t30 * t58 - t8 * t36) + g(3) * (-t18 * t36 - t31 * t33);
t19 = t21 * t29;
t40 = -g(1) * t12 + g(2) * t7 + g(3) * t19;
t38 = -g(1) * t43 - g(3) * t56;
t22 = pkin(2) * t46;
t17 = -g(3) * t31 + (-g(1) * t28 + g(2) * t30) * t29;
t14 = -t18 * t33 + t31 * t36;
t5 = t30 * t57 + t8 * t33;
t3 = -t11 * t33 + t28 * t57;
t1 = g(1) * t11 + g(2) * t8 + g(3) * t18;
t2 = [-g(3), 0, 0, 0, -g(3), 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(2) * (t46 - t59) + t38, -g(1) * (t28 * t54 - t30 * t37) - g(2) * (-t28 * t37 - t30 * t54) + g(3) * t29 * t34, -g(2) * t22 + (g(2) * t59 + t38) * pkin(2), t1, t40, -g(1) * (t43 * pkin(2) + t11 * pkin(3) + qJ(4) * t12) - g(2) * (-pkin(2) * t59 + t8 * pkin(3) - t7 * qJ(4) + t22) - g(3) * (pkin(2) * t56 + t18 * pkin(3) - t19 * qJ(4)) 0, 0, 0, 0, 0, t40 * t33, t40 * t36, 0, 0, 0, 0, 0, -g(1) * (t11 * t32 + t12 * t51) - g(2) * (t8 * t32 - t7 * t51) - g(3) * (t18 * t32 - t19 * t51) -g(1) * (t11 * t35 - t12 * t52) - g(2) * (t8 * t35 + t7 * t52) - g(3) * (t18 * t35 + t19 * t52); 0, 0, 0, 0, t17, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, g(1) * t3 - g(2) * t5 + g(3) * t14, 0, 0, 0, 0, 0, -t41 * t35, t41 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t35 - t3 * t32) - g(2) * (t5 * t32 - t35 * t7) - g(3) * (-t14 * t32 - t19 * t35) -g(1) * (-t12 * t32 - t3 * t35) - g(2) * (t32 * t7 + t5 * t35) - g(3) * (-t14 * t35 + t19 * t32);];
taug_reg  = t2;
