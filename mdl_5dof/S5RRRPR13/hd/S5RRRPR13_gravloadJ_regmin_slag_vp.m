% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:33
% EndTime: 2019-12-31 21:46:34
% DurationCPUTime: 0.35s
% Computational Cost: add. (225->78), mult. (600->127), div. (0->0), fcn. (743->10), ass. (0->46)
t26 = sin(qJ(2));
t27 = sin(qJ(1));
t30 = cos(qJ(2));
t43 = cos(pkin(5));
t53 = cos(qJ(1));
t36 = t43 * t53;
t14 = t27 * t26 - t30 * t36;
t24 = sin(qJ(5));
t28 = cos(qJ(5));
t15 = t26 * t36 + t27 * t30;
t25 = sin(qJ(3));
t29 = cos(qJ(3));
t23 = sin(pkin(5));
t42 = t23 * t53;
t6 = t15 * t25 + t29 * t42;
t58 = t14 * t28 + t6 * t24;
t57 = -t14 * t24 + t6 * t28;
t41 = t27 * t43;
t16 = t53 * t26 + t30 * t41;
t54 = g(3) * t23;
t32 = -g(1) * t16 - g(2) * t14 + t30 * t54;
t50 = t23 * t26;
t49 = t23 * t27;
t48 = t23 * t29;
t47 = t24 * t25;
t46 = t24 * t30;
t45 = t25 * t28;
t44 = t28 * t30;
t7 = t15 * t29 - t25 * t42;
t17 = -t26 * t41 + t53 * t30;
t10 = t17 * t25 - t27 * t48;
t40 = -g(1) * t6 + g(2) * t10;
t11 = t17 * t29 + t25 * t49;
t39 = -g(1) * t7 + g(2) * t11;
t38 = g(1) * t14 - g(2) * t16;
t37 = -g(1) * t17 - g(2) * t15;
t12 = t25 * t50 - t43 * t29;
t34 = g(1) * t10 + g(2) * t6 + g(3) * t12;
t13 = t43 * t25 + t26 * t48;
t33 = g(1) * t11 + g(2) * t7 + g(3) * t13;
t31 = g(3) * t50 - t37;
t5 = t32 * t29;
t4 = t32 * t25;
t3 = t10 * t24 + t16 * t28;
t2 = t10 * t28 - t16 * t24;
t1 = [0, g(1) * t27 - g(2) * t53, g(1) * t53 + g(2) * t27, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t17, -t38, 0, 0, 0, 0, 0, -t39, t40, t38, t39, -t40, -g(1) * (-t27 * pkin(1) - t15 * pkin(2) - pkin(3) * t7 + pkin(7) * t42 - t14 * pkin(8) - qJ(4) * t6) - g(2) * (t53 * pkin(1) + t17 * pkin(2) + t11 * pkin(3) + pkin(7) * t49 + t16 * pkin(8) + t10 * qJ(4)), 0, 0, 0, 0, 0, g(1) * t58 - g(2) * t3, g(1) * t57 - g(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, -t32, t31, 0, 0, 0, 0, 0, -t5, t4, -t31, t5, -t4, (-t26 * t54 + t37) * pkin(8) - t32 * (pkin(3) * t29 + qJ(4) * t25 + pkin(2)), 0, 0, 0, 0, 0, -g(1) * (-t16 * t47 + t17 * t28) - g(2) * (-t14 * t47 + t15 * t28) - (t25 * t46 + t26 * t28) * t54, -g(1) * (-t16 * t45 - t17 * t24) - g(2) * (-t14 * t45 - t15 * t24) - (-t24 * t26 + t25 * t44) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t33, 0, -t34, -t33, -g(1) * (-t10 * pkin(3) + t11 * qJ(4)) - g(2) * (-t6 * pkin(3) + t7 * qJ(4)) - g(3) * (-t12 * pkin(3) + t13 * qJ(4)), 0, 0, 0, 0, 0, -t33 * t24, -t33 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t57 - g(3) * (t12 * t28 + t23 * t46), g(1) * t3 + g(2) * t58 - g(3) * (-t12 * t24 + t23 * t44);];
taug_reg = t1;
