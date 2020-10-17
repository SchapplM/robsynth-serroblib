% Calculate inertial parameters regressor of gravitation load for
% S5RRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:28
% EndTime: 2019-12-31 20:07:28
% DurationCPUTime: 0.30s
% Computational Cost: add. (237->71), mult. (372->98), div. (0->0), fcn. (382->8), ass. (0->43)
t29 = cos(pkin(8));
t20 = t29 * pkin(3) + pkin(2);
t33 = cos(qJ(2));
t13 = t33 * t20;
t30 = -pkin(7) - qJ(3);
t31 = sin(qJ(2));
t58 = -t31 * t30 + t13;
t32 = sin(qJ(1));
t34 = cos(qJ(1));
t14 = g(1) * t34 + g(2) * t32;
t9 = -g(3) * t33 + t14 * t31;
t57 = g(1) * t32;
t54 = g(3) * t31;
t28 = sin(pkin(8));
t51 = t32 * t28;
t50 = t32 * t33;
t27 = pkin(8) + qJ(4);
t21 = sin(t27);
t49 = t34 * t21;
t22 = cos(t27);
t48 = t34 * t22;
t47 = t34 * t28;
t46 = t34 * t29;
t45 = t34 * pkin(1) + t32 * pkin(6);
t5 = t21 * t50 + t48;
t7 = -t32 * t22 + t33 * t49;
t44 = g(1) * t5 - g(2) * t7;
t43 = -g(2) * t34 + t57;
t42 = t33 * pkin(2) + t31 * qJ(3);
t40 = pkin(4) * t22 + qJ(5) * t21;
t37 = pkin(3) * t51 + t58 * t34 + t45;
t1 = g(1) * t7 + g(2) * t5 + t21 * t54;
t6 = t22 * t50 - t49;
t8 = t32 * t21 + t33 * t48;
t36 = g(1) * t8 + g(2) * t6 + t22 * t54;
t24 = t34 * pkin(6);
t35 = pkin(3) * t47 + t24 + (-pkin(1) - t58) * t32;
t11 = t43 * t31;
t10 = t14 * t33 + t54;
t4 = t9 * t22;
t3 = t9 * t21;
t2 = g(1) * t6 - g(2) * t8;
t12 = [0, 0, 0, 0, 0, 0, t43, t14, 0, 0, 0, 0, 0, 0, 0, 0, t43 * t33, -t11, -t14, -g(1) * (-t32 * pkin(1) + t24) - g(2) * t45, 0, 0, 0, 0, 0, 0, -g(1) * (-t29 * t50 + t47) - g(2) * (t33 * t46 + t51), -g(1) * (t28 * t50 + t46) - g(2) * (t32 * t29 - t33 * t47), t11, -g(1) * t24 - g(2) * (t42 * t34 + t45) - (-pkin(1) - t42) * t57, 0, 0, 0, 0, 0, 0, t2, -t44, t11, -g(1) * t35 - g(2) * t37, 0, 0, 0, 0, 0, 0, t2, t11, t44, -g(1) * (-t6 * pkin(4) - t5 * qJ(5) + t35) - g(2) * (t8 * pkin(4) + t7 * qJ(5) + t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t29, -t9 * t28, -t10, -g(3) * t42 + t14 * (pkin(2) * t31 - qJ(3) * t33), 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(3) * t58 + t14 * (t20 * t31 + t30 * t33), 0, 0, 0, 0, 0, 0, t4, -t10, t3, -g(3) * t13 + (-g(3) * t40 + t14 * t30) * t33 + (g(3) * t30 + t14 * (t20 + t40)) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t36, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t36, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (-t5 * pkin(4) + t6 * qJ(5)) - (-pkin(4) * t21 + qJ(5) * t22) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t12;
