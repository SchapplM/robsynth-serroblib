% Calculate minimal parameter regressor of gravitation load for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% taug_reg [6x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:10:55
% EndTime: 2019-05-04 22:10:56
% DurationCPUTime: 0.33s
% Computational Cost: add. (257->77), mult. (592->141), div. (0->0), fcn. (750->14), ass. (0->49)
t28 = sin(pkin(11));
t36 = sin(qJ(2));
t39 = cos(qJ(2));
t49 = cos(pkin(11));
t45 = -t36 * t28 + t39 * t49;
t27 = qJ(4) + pkin(12);
t26 = cos(t27);
t34 = sin(qJ(6));
t60 = t26 * t34;
t37 = cos(qJ(6));
t59 = t26 * t37;
t29 = sin(pkin(10));
t30 = sin(pkin(6));
t58 = t29 * t30;
t57 = t29 * t36;
t31 = cos(pkin(10));
t56 = t30 * t31;
t35 = sin(qJ(4));
t55 = t30 * t35;
t38 = cos(qJ(4));
t54 = t30 * t38;
t53 = t30 * t39;
t32 = cos(pkin(6));
t52 = t32 * t36;
t51 = t32 * t39;
t48 = t31 * t51;
t20 = -t39 * t28 - t36 * t49;
t18 = t20 * t32;
t7 = -t31 * t18 + t29 * t45;
t8 = -t29 * t18 - t31 * t45;
t46 = -t29 * t51 - t31 * t36;
t17 = t20 * t30;
t25 = sin(t27);
t44 = g(1) * (t25 * t8 + t26 * t58) + g(2) * (-t7 * t25 - t26 * t56) + g(3) * (t17 * t25 + t32 * t26);
t16 = t45 * t30;
t42 = t45 * t32;
t6 = t29 * t20 + t31 * t42;
t9 = t31 * t20 - t29 * t42;
t43 = g(1) * t9 + g(2) * t6 + g(3) * t16;
t41 = -g(1) * t46 - g(3) * t53;
t40 = -g(1) * (t29 * t54 + t35 * t8) - g(2) * (-t31 * t54 - t7 * t35) - g(3) * (t17 * t35 + t32 * t38);
t33 = -qJ(5) - pkin(8);
t24 = t38 * pkin(4) + pkin(3);
t21 = pkin(2) * t48;
t15 = -g(3) * t32 + (-g(1) * t29 + g(2) * t31) * t30;
t12 = -t17 * t26 + t32 * t25;
t4 = t25 * t58 - t26 * t8;
t2 = -t25 * t56 + t7 * t26;
t1 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(2) * (t48 - t57) + t41, -g(1) * (t29 * t52 - t31 * t39) - g(2) * (-t29 * t39 - t31 * t52) + g(3) * t30 * t36, -g(2) * t21 + (g(2) * t57 + t41) * pkin(2), 0, 0, 0, 0, 0, -t43 * t38, t43 * t35, g(1) * t8 - g(2) * t7 + g(3) * t17, -g(1) * (t46 * pkin(2) + t9 * t24 + t8 * t33) - g(2) * (-pkin(2) * t57 + t6 * t24 - t33 * t7 + t21) - g(3) * (pkin(2) * t53 + t16 * t24 + t17 * t33) 0, 0, 0, 0, 0, -g(1) * (-t8 * t34 + t9 * t59) - g(2) * (t34 * t7 + t6 * t59) - g(3) * (t16 * t59 - t17 * t34) -g(1) * (-t8 * t37 - t9 * t60) - g(2) * (t37 * t7 - t6 * t60) - g(3) * (-t16 * t60 - t17 * t37); 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -g(1) * (-t29 * t55 + t38 * t8) - g(2) * (t31 * t55 - t7 * t38) - g(3) * (t17 * t38 - t32 * t35) 0, t40 * pkin(4), 0, 0, 0, 0, 0, -t44 * t37, t44 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t34 - t9 * t37) - g(2) * (-t2 * t34 - t6 * t37) - g(3) * (-t12 * t34 - t16 * t37) -g(1) * (t9 * t34 - t4 * t37) - g(2) * (-t2 * t37 + t6 * t34) - g(3) * (-t12 * t37 + t16 * t34);];
taug_reg  = t1;
