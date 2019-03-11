% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t15 = g(1) * t35 + g(2) * t32;
t28 = qJ(4) + qJ(5);
t21 = sin(t28);
t27 = qJ(2) + pkin(10);
t20 = cos(t27);
t22 = cos(t28);
t45 = t35 * t22;
t50 = t32 * t21;
t3 = t20 * t50 + t45;
t46 = t35 * t21;
t49 = t32 * t22;
t5 = -t20 * t46 + t49;
t19 = sin(t27);
t53 = g(3) * t19;
t1 = -g(1) * t5 + g(2) * t3 + t21 * t53;
t38 = -g(3) * t20 + t15 * t19;
t30 = sin(qJ(4));
t12 = t30 * pkin(4) + pkin(5) * t21;
t51 = t12 * t20;
t48 = t32 * t30;
t33 = cos(qJ(4));
t47 = t32 * t33;
t44 = t35 * t30;
t43 = t35 * t33;
t29 = -qJ(3) - pkin(7);
t42 = t12 - t29;
t13 = t33 * pkin(4) + pkin(5) * t22;
t14 = g(1) * t32 - g(2) * t35;
t11 = pkin(3) + t13;
t26 = -qJ(6) - pkin(9) - pkin(8);
t40 = t20 * t11 - t19 * t26;
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t36 = -g(3) * t34 + t15 * t31;
t24 = t34 * pkin(2);
t18 = t24 + pkin(1);
t16 = t35 * t18;
t10 = t20 * t43 + t48;
t9 = -t20 * t44 + t47;
t8 = -t20 * t47 + t44;
t7 = t20 * t48 + t43;
t6 = t20 * t45 + t50;
t4 = -t20 * t49 + t46;
t2 = g(1) * t6 - g(2) * t4 + t22 * t53;
t17 = [0, t14, t15, 0, 0, 0, 0, 0, t14 * t34, -t14 * t31, -t15, -g(1) * (-t32 * t18 - t35 * t29) - g(2) * (-t32 * t29 + t16) 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t14 * t19, -g(2) * t16 + (-g(1) * t42 - g(2) * t40) * t35 + (-g(1) * (-t18 - t40) - g(2) * t42) * t32; 0, 0, 0, 0, 0, 0, 0, 0, t36, g(3) * t31 + t15 * t34, 0, t36 * pkin(2), 0, 0, 0, 0, 0, t38 * t33, -t38 * t30, 0, 0, 0, 0, 0, t38 * t22, -t38 * t21, -t15 * t20 - t53, -g(3) * (t24 + t40) + t15 * (pkin(2) * t31 + t11 * t19 + t20 * t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t30 * t53, g(1) * t10 - g(2) * t8 + t33 * t53, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t32 * t13 - t35 * t51) - g(2) * (-t35 * t13 - t32 * t51) + t12 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38;];
taug_reg  = t17;
