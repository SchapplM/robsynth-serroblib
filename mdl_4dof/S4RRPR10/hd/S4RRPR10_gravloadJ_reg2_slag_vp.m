% Calculate inertial parameters regressor of gravitation load for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t21 = sin(qJ(2));
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t9 = g(1) * t25 + g(2) * t22;
t46 = t9 * t21;
t14 = t21 * qJ(3);
t24 = cos(qJ(2));
t32 = t24 * pkin(2) + t14;
t2 = g(3) * t21 + t9 * t24;
t44 = pkin(2) * t21;
t43 = g(1) * t22;
t39 = g(3) * t24;
t38 = t24 * pkin(6);
t20 = sin(qJ(4));
t37 = t22 * t20;
t23 = cos(qJ(4));
t36 = t22 * t23;
t35 = t24 * t25;
t34 = t25 * t20;
t33 = t25 * t23;
t31 = t25 * pkin(1) + t22 * pkin(5);
t30 = qJ(3) * t24;
t29 = pkin(2) * t35 + t25 * t14 + t31;
t28 = -g(2) * t25 + t43;
t27 = -pkin(1) - t32;
t17 = t25 * pkin(5);
t12 = t25 * t30;
t10 = t22 * t30;
t8 = t28 * t24;
t7 = t28 * t21;
t6 = -t21 * t37 + t33;
t5 = t21 * t36 + t34;
t4 = t21 * t34 + t36;
t3 = t21 * t33 - t37;
t1 = -t39 + t46;
t11 = [0, 0, 0, 0, 0, 0, t28, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t9, -g(1) * (-t22 * pkin(1) + t17) - g(2) * t31, 0, 0, 0, 0, 0, 0, -t9, -t8, t7, -g(1) * t17 - g(2) * t29 - t27 * t43, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t8, -g(1) * (t25 * pkin(3) + t17) - g(2) * (pkin(6) * t35 + t29) + (-g(1) * (t27 - t38) - g(2) * pkin(3)) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (-t25 * t44 + t12) - g(2) * (-t22 * t44 + t10) - g(3) * t32, 0, 0, 0, 0, 0, 0, -t2 * t20, -t2 * t23, t1, -g(1) * t12 - g(2) * t10 - g(3) * (t32 + t38) + (pkin(2) + pkin(6)) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t23 * t39, g(1) * t4 - g(2) * t6 - t20 * t39, 0, 0;];
taug_reg = t11;
