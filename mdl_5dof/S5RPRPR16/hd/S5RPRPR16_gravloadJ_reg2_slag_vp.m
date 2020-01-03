% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR16_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t24 = sin(qJ(3));
t27 = cos(qJ(3));
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t44 = g(2) * t28;
t9 = g(1) * t25 - t44;
t1 = g(3) * t27 + t9 * t24;
t46 = -pkin(1) - pkin(6);
t45 = -pkin(3) - pkin(7);
t43 = g(3) * t24;
t41 = t24 * pkin(3);
t40 = t24 * t28;
t39 = t25 * t27;
t23 = sin(qJ(5));
t38 = t28 * t23;
t26 = cos(qJ(5));
t37 = t28 * t26;
t34 = qJ(4) * t24;
t36 = pkin(3) * t39 + t25 * t34;
t35 = t28 * pkin(1) + t25 * qJ(2);
t17 = t27 * qJ(4);
t33 = t28 * pkin(6) + t35;
t32 = t46 * t25;
t31 = t25 * t41 + t33;
t10 = g(1) * t28 + g(2) * t25;
t18 = t28 * qJ(2);
t30 = pkin(3) * t40 - t28 * t17 + t18;
t8 = t10 * t27;
t7 = t10 * t24;
t6 = -t23 * t39 + t37;
t5 = -t26 * t39 - t38;
t4 = -t25 * t26 - t27 * t38;
t3 = t25 * t23 - t27 * t37;
t2 = g(1) * t39 - t27 * t44 - t43;
t11 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, -g(1) * (-t25 * pkin(1) + t18) - g(2) * t35, 0, 0, 0, 0, 0, 0, -t7, -t8, t9, -g(1) * (t18 + t32) - g(2) * t33, 0, 0, 0, 0, 0, 0, t9, t7, t8, -g(1) * (t32 + t30) - g(2) * (-t25 * t17 + t31), 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, -t7, -g(1) * (pkin(7) * t40 + t30) - g(2) * (t28 * pkin(4) + t31) + (-g(1) * (-pkin(4) + t46) - g(2) * (t24 * pkin(7) - t17)) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -g(1) * t36 - g(3) * (t17 - t41) - (-pkin(3) * t27 - t34) * t44, 0, 0, 0, 0, 0, 0, -t1 * t23, -t1 * t26, -t2, -g(1) * (pkin(7) * t39 + t36) - g(3) * (t45 * t24 + t17) - (t45 * t27 - t34) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 - t26 * t43, g(1) * t6 - g(2) * t4 + t23 * t43, 0, 0;];
taug_reg = t11;
