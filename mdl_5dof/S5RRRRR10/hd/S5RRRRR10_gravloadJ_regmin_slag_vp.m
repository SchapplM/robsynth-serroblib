% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:35:58
% EndTime: 2019-12-31 22:35:59
% DurationCPUTime: 0.33s
% Computational Cost: add. (270->66), mult. (502->123), div. (0->0), fcn. (634->12), ass. (0->47)
t29 = sin(qJ(2));
t30 = sin(qJ(1));
t33 = cos(qJ(2));
t34 = cos(qJ(1));
t42 = cos(pkin(5));
t39 = t34 * t42;
t16 = t30 * t29 - t33 * t39;
t27 = sin(qJ(5));
t31 = cos(qJ(5));
t17 = t29 * t39 + t30 * t33;
t25 = qJ(3) + qJ(4);
t23 = sin(t25);
t24 = cos(t25);
t26 = sin(pkin(5));
t45 = t26 * t34;
t8 = t17 * t24 - t23 * t45;
t55 = -t16 * t31 + t8 * t27;
t54 = t16 * t27 + t8 * t31;
t53 = g(3) * t26;
t50 = t24 * t27;
t49 = t24 * t31;
t48 = t26 * t29;
t47 = t26 * t30;
t32 = cos(qJ(3));
t46 = t26 * t32;
t44 = t27 * t33;
t43 = t31 * t33;
t28 = sin(qJ(3));
t41 = t17 * t32 - t28 * t45;
t40 = t30 * t42;
t38 = t17 * t23 + t24 * t45;
t37 = t17 * t28 + t32 * t45;
t19 = -t29 * t40 + t34 * t33;
t10 = -t19 * t23 + t24 * t47;
t36 = g(1) * t10 - g(2) * t38 + g(3) * (-t23 * t48 + t42 * t24);
t18 = t34 * t29 + t33 * t40;
t35 = -g(1) * t18 - g(2) * t16 + t33 * t53;
t15 = t42 * t23 + t24 * t48;
t13 = t19 * t32 + t28 * t47;
t12 = -t19 * t28 + t30 * t46;
t11 = t19 * t24 + t23 * t47;
t6 = t11 * t31 + t18 * t27;
t5 = -t11 * t27 + t18 * t31;
t4 = g(1) * t11 + g(2) * t8 + g(3) * t15;
t2 = t36 * t31;
t1 = t36 * t27;
t3 = [0, g(1) * t30 - g(2) * t34, g(1) * t34 + g(2) * t30, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t19, -g(1) * t16 + g(2) * t18, 0, 0, 0, 0, 0, g(1) * t41 - g(2) * t13, -g(1) * t37 - g(2) * t12, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t11, -g(1) * t38 - g(2) * t10, 0, 0, 0, 0, 0, g(1) * t54 - g(2) * t6, -g(1) * t55 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, -t35, g(1) * t19 + g(2) * t17 + g(3) * t48, 0, 0, 0, 0, 0, -t35 * t32, t35 * t28, 0, 0, 0, 0, 0, -t35 * t24, t35 * t23, 0, 0, 0, 0, 0, -g(1) * (-t18 * t49 + t19 * t27) - g(2) * (-t16 * t49 + t17 * t27) - (t24 * t43 + t27 * t29) * t53, -g(1) * (t18 * t50 + t19 * t31) - g(2) * (t16 * t50 + t17 * t31) - (-t24 * t44 + t29 * t31) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 + g(2) * t37 - g(3) * (-t28 * t48 + t42 * t32), g(1) * t13 + g(2) * t41 - g(3) * (-t42 * t28 - t29 * t46), 0, 0, 0, 0, 0, -t36, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t55 - g(3) * (-t15 * t27 - t26 * t43), g(1) * t6 + g(2) * t54 - g(3) * (-t15 * t31 + t26 * t44);];
taug_reg = t3;
