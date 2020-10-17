% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:22
% EndTime: 2020-01-03 11:36:23
% DurationCPUTime: 0.22s
% Computational Cost: add. (236->46), mult. (158->58), div. (0->0), fcn. (158->10), ass. (0->36)
t29 = sin(pkin(9));
t30 = cos(pkin(9));
t49 = pkin(4) * t30 + pkin(7) * t29;
t28 = qJ(1) + pkin(8);
t25 = qJ(3) + t28;
t21 = sin(t25);
t48 = t49 * t21;
t22 = cos(t25);
t47 = t49 * t22;
t44 = g(1) * t29;
t31 = sin(qJ(5));
t43 = t30 * t31;
t33 = cos(qJ(5));
t42 = t30 * t33;
t41 = t22 * pkin(3) + t21 * qJ(4);
t23 = sin(t28);
t32 = sin(qJ(1));
t40 = t32 * pkin(1) + pkin(2) * t23;
t24 = cos(t28);
t34 = cos(qJ(1));
t39 = t34 * pkin(1) + pkin(2) * t24;
t38 = t21 * pkin(3) - t22 * qJ(4);
t37 = t39 + t41;
t10 = g(2) * t22 + g(3) * t21;
t36 = -g(2) * t34 - g(3) * t32;
t35 = t38 + t40;
t9 = g(2) * t21 - g(3) * t22;
t8 = t10 * t30;
t7 = t10 * t29;
t6 = t21 * t31 + t22 * t42;
t5 = -t21 * t33 + t22 * t43;
t4 = t21 * t42 - t22 * t31;
t3 = -t21 * t43 - t22 * t33;
t2 = -g(2) * t6 - g(3) * t4;
t1 = g(2) * t5 - g(3) * t3;
t11 = [0, 0, 0, 0, 0, 0, t36, g(2) * t32 - g(3) * t34, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t24 - g(3) * t23, g(2) * t23 - g(3) * t24, 0, t36 * pkin(1), 0, 0, 0, 0, 0, 0, -t10, t9, 0, -g(2) * t39 - g(3) * t40, 0, 0, 0, 0, 0, 0, -t8, t7, -t9, -g(2) * t37 - g(3) * t35, 0, 0, 0, 0, 0, 0, t2, t1, -t7, -g(2) * (t37 + t47) - g(3) * (t35 + t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, -t9, -g(2) * t41 - g(3) * t38, 0, 0, 0, 0, 0, 0, t2, t1, -t7, -g(2) * (t41 + t47) - g(3) * (t38 + t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t3 - g(3) * t5 + t31 * t44, g(2) * t4 - g(3) * t6 + t33 * t44, 0, 0;];
taug_reg = t11;
