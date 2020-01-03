% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t25 = qJ(3) + qJ(4);
t19 = sin(t25);
t20 = cos(t25);
t38 = -t19 * pkin(4) + t20 * pkin(8);
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t51 = g(2) * t31;
t10 = g(1) * t28 - t51;
t54 = -g(3) * t19 + t10 * t20;
t30 = cos(qJ(3));
t53 = pkin(3) * t30;
t52 = pkin(4) * t20;
t49 = g(3) * t20;
t27 = sin(qJ(3));
t47 = t27 * pkin(3);
t46 = t19 * t28;
t26 = sin(qJ(5));
t45 = t28 * t26;
t29 = cos(qJ(5));
t44 = t28 * t29;
t43 = t31 * t26;
t42 = t31 * t29;
t41 = pkin(8) * t46 + t28 * t52;
t40 = t31 * pkin(1) + t28 * qJ(2);
t22 = t31 * qJ(2);
t39 = -t28 * pkin(1) + t22;
t37 = -pkin(8) * t19 - t52;
t11 = g(1) * t31 + g(2) * t28;
t32 = -pkin(7) - pkin(6);
t35 = t28 * t32 + t31 * t47 + t39;
t34 = t28 * t47 - t31 * t32 + t40;
t33 = g(3) * t27 - t10 * t30;
t9 = t19 * t42 - t45;
t8 = t19 * t43 + t44;
t7 = t19 * t44 + t43;
t6 = -t19 * t45 + t42;
t5 = t11 * t20;
t3 = g(1) * t46 - t19 * t51 + t49;
t2 = t54 * t29;
t1 = t54 * t26;
t4 = [0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, -g(1) * t39 - g(2) * t40, 0, 0, 0, 0, 0, 0, -t11 * t27, -t11 * t30, t10, -g(1) * (t22 + (-pkin(1) - pkin(6)) * t28) - g(2) * (t31 * pkin(6) + t40), 0, 0, 0, 0, 0, 0, -t11 * t19, -t5, t10, -g(1) * t35 - g(2) * t34, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7, g(1) * t8 - g(2) * t6, t5, -g(1) * (-t31 * t38 + t35) - g(2) * (-t28 * t38 + t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, g(3) * t30 + t10 * t27, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t3, 0, t33 * pkin(3), 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * (t28 * t53 + t41) - g(3) * (t38 - t47) - (t37 - t53) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * t41 - g(3) * t38 - t37 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8 + t26 * t49, g(1) * t7 - g(2) * t9 + t29 * t49, 0, 0;];
taug_reg = t4;
