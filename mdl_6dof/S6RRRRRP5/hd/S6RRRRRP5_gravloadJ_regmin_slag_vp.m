% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 05:03:21
% EndTime: 2019-05-08 05:03:22
% DurationCPUTime: 0.29s
% Computational Cost: add. (349->64), mult. (362->105), div. (0->0), fcn. (384->10), ass. (0->52)
t37 = sin(qJ(1));
t40 = cos(qJ(1));
t46 = g(1) * t40 + g(2) * t37;
t34 = qJ(3) + qJ(4);
t32 = qJ(5) + t34;
t27 = sin(t32);
t28 = cos(t32);
t52 = t40 * t28;
t39 = cos(qJ(2));
t55 = t37 * t39;
t5 = t27 * t55 + t52;
t36 = sin(qJ(2));
t58 = g(3) * t36;
t53 = t40 * t27;
t7 = t28 * t37 - t39 * t53;
t1 = -g(1) * t7 + g(2) * t5 + t27 * t58;
t14 = -g(3) * t39 + t36 * t46;
t29 = sin(t34);
t23 = -pkin(4) * t29 - pkin(5) * t27;
t35 = sin(qJ(3));
t16 = pkin(3) * t35 - t23;
t56 = pkin(7) + t16;
t54 = t39 * t40;
t51 = t40 * t29;
t30 = cos(t34);
t50 = t40 * t30;
t49 = t40 * t35;
t38 = cos(qJ(3));
t48 = t40 * t38;
t24 = pkin(4) * t30 + pkin(5) * t28;
t17 = pkin(3) * t38 + t24;
t45 = g(1) * t37 - g(2) * t40;
t31 = -qJ(6) - pkin(10) - pkin(9) - pkin(8);
t9 = pkin(2) + t17;
t43 = -t31 * t36 + t39 * t9;
t42 = -pkin(1) - t43;
t22 = t45 * t36;
t21 = t35 * t37 + t39 * t48;
t20 = t37 * t38 - t39 * t49;
t19 = -t38 * t55 + t49;
t18 = t35 * t55 + t48;
t15 = t39 * t46 + t58;
t13 = t29 * t37 + t39 * t50;
t12 = t30 * t37 - t39 * t51;
t11 = -t30 * t55 + t51;
t10 = t29 * t55 + t50;
t8 = t27 * t37 + t39 * t52;
t6 = -t28 * t55 + t53;
t4 = g(1) * t13 - g(2) * t11 + t30 * t58;
t3 = -g(1) * t12 + g(2) * t10 + t29 * t58;
t2 = g(1) * t8 - g(2) * t6 + t28 * t58;
t25 = [0, t45, t46, 0, 0, 0, 0, 0, t45 * t39, -t22, 0, 0, 0, 0, 0, -g(1) * t19 - g(2) * t21, -g(1) * t18 - g(2) * t20, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13, -g(1) * t10 - g(2) * t12, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, t22 (-g(1) * t56 + g(2) * t42) * t40 + (-g(1) * t42 - g(2) * t56) * t37; 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t14 * t38, -t14 * t35, 0, 0, 0, 0, 0, t14 * t30, -t14 * t29, 0, 0, 0, 0, 0, t14 * t28, -t14 * t27, -t15, -g(3) * t43 + t46 * (t31 * t39 + t36 * t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t20 + g(2) * t18 + t35 * t58, g(1) * t21 - g(2) * t19 + t38 * t58, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t16 * t54 + t17 * t37) - g(2) * (-t16 * t55 - t17 * t40) + t16 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t23 * t54 + t24 * t37) - g(2) * (t23 * t55 - t24 * t40) - t23 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14;];
taug_reg  = t25;
