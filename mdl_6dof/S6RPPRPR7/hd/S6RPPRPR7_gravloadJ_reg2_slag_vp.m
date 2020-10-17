% Calculate inertial parameters regressor of gravitation load for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:32:14
% EndTime: 2019-05-05 14:32:15
% DurationCPUTime: 0.33s
% Computational Cost: add. (222->71), mult. (269->89), div. (0->0), fcn. (265->10), ass. (0->44)
t30 = -pkin(7) - qJ(3);
t31 = sin(qJ(1));
t32 = cos(qJ(1));
t26 = sin(pkin(9));
t56 = pkin(3) * t26;
t58 = t31 * t30 + t32 * t56;
t57 = -g(1) * t31 + g(2) * t32;
t24 = pkin(9) + qJ(4);
t16 = sin(t24);
t18 = cos(t24);
t2 = -g(3) * t16 - t18 * t57;
t25 = sin(pkin(10));
t55 = pkin(5) * t25;
t51 = g(3) * t18;
t23 = pkin(10) + qJ(6);
t15 = sin(t23);
t50 = t31 * t15;
t17 = cos(t23);
t49 = t31 * t17;
t48 = t31 * t25;
t27 = cos(pkin(10));
t47 = t31 * t27;
t46 = t32 * t15;
t45 = t32 * t17;
t44 = t32 * t25;
t43 = t32 * t27;
t42 = t32 * pkin(1) + t31 * qJ(2);
t41 = t31 * t56 + t42;
t20 = t32 * qJ(2);
t40 = -t31 * pkin(1) + t20;
t9 = g(1) * t32 + g(2) * t31;
t38 = t16 * pkin(4) - t18 * qJ(5);
t13 = t27 * pkin(5) + pkin(4);
t29 = -pkin(8) - qJ(5);
t36 = t16 * t13 + t18 * t29;
t35 = t40 + t58;
t34 = -t32 * t30 + t41;
t7 = t9 * t18;
t6 = t16 * t45 - t50;
t5 = t16 * t46 + t49;
t4 = t16 * t49 + t46;
t3 = -t16 * t50 + t45;
t1 = -t57 * t16 + t51;
t8 = [0, 0, 0, 0, 0, 0, -t57, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t9, -g(1) * t40 - g(2) * t42, 0, 0, 0, 0, 0, 0, -t9 * t26, -t9 * cos(pkin(9)) -t57, -g(1) * (t20 + (-pkin(1) - qJ(3)) * t31) - g(2) * (t32 * qJ(3) + t42) 0, 0, 0, 0, 0, 0, -t9 * t16, -t7, -t57, -g(1) * t35 - g(2) * t34, 0, 0, 0, 0, 0, 0, -g(1) * (t16 * t43 - t48) - g(2) * (t16 * t47 + t44) -g(1) * (-t16 * t44 - t47) - g(2) * (-t16 * t48 + t43) t7, -g(1) * (t32 * t38 + t35) - g(2) * (t31 * t38 + t34) 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t7, -g(1) * (t20 + t58) - g(2) * t41 + (-g(1) * t36 - g(2) * (-t30 + t55)) * t32 + (-g(1) * (-pkin(1) - t55) - g(2) * t36) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t27, t2 * t25, -t1, g(3) * t38 + t57 * (pkin(4) * t18 + qJ(5) * t16) 0, 0, 0, 0, 0, 0, -t2 * t17, t2 * t15, -t1, g(3) * t36 + t57 * (t13 * t18 - t16 * t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t15 * t51, g(1) * t4 - g(2) * t6 + t17 * t51, 0, 0;];
taug_reg  = t8;
