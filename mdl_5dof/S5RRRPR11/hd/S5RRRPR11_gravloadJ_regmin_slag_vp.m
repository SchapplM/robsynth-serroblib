% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t20 = sin(qJ(2));
t24 = cos(qJ(2));
t21 = sin(qJ(1));
t25 = cos(qJ(1));
t34 = g(1) * t25 + g(2) * t21;
t54 = g(3) * t24 - t34 * t20;
t23 = cos(qJ(3));
t19 = sin(qJ(3));
t38 = t25 * t19;
t13 = -t21 * t23 + t24 * t38;
t37 = t25 * t23;
t14 = t21 * t19 + t24 * t37;
t18 = sin(qJ(5));
t22 = cos(qJ(5));
t2 = t13 * t22 - t14 * t18;
t31 = t18 * t23 - t19 * t22;
t39 = t21 * t24;
t11 = t19 * t39 + t37;
t12 = t23 * t39 - t38;
t36 = -t11 * t22 + t12 * t18;
t42 = g(3) * t20;
t53 = -g(1) * t2 + g(2) * t36 + t31 * t42;
t3 = t13 * t18 + t14 * t22;
t30 = t18 * t19 + t22 * t23;
t32 = t11 * t18 + t12 * t22;
t49 = g(1) * t3 + g(2) * t32 + t30 * t42;
t35 = g(1) * t11 - g(2) * t13;
t33 = g(1) * t21 - g(2) * t25;
t29 = t24 * pkin(2) + t20 * pkin(7) + pkin(1);
t28 = pkin(3) * t23 + qJ(4) * t19 + pkin(2);
t1 = g(1) * t13 + g(2) * t11 + t19 * t42;
t27 = g(1) * t14 + g(2) * t12 + t23 * t42;
t15 = t33 * t20;
t7 = t34 * t24 + t42;
t6 = t54 * t23;
t5 = t54 * t19;
t4 = g(1) * t12 - g(2) * t14;
t8 = [0, t33, t34, 0, 0, 0, 0, 0, t33 * t24, -t15, 0, 0, 0, 0, 0, t4, -t35, t4, t15, t35, -g(1) * (-t12 * pkin(3) - t11 * qJ(4)) - g(2) * (t14 * pkin(3) + t13 * qJ(4)) + (-g(1) * pkin(6) - g(2) * t29) * t25 + (-g(2) * pkin(6) + g(1) * t29) * t21, 0, 0, 0, 0, 0, g(1) * t32 - g(2) * t3, -g(1) * t36 - g(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, -t54, t7, 0, 0, 0, 0, 0, -t6, t5, -t6, -t7, -t5, (-t34 * pkin(7) - g(3) * t28) * t24 + (-g(3) * pkin(7) + t34 * t28) * t20, 0, 0, 0, 0, 0, -t54 * t30, t54 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t27, t1, 0, -t27, -g(1) * (-t13 * pkin(3) + t14 * qJ(4)) - g(2) * (-t11 * pkin(3) + t12 * qJ(4)) - (-pkin(3) * t19 + qJ(4) * t23) * t42, 0, 0, 0, 0, 0, -t53, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t49;];
taug_reg = t8;
