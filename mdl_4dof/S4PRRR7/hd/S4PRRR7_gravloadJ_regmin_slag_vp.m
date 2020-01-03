% Calculate minimal parameter regressor of gravitation load for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [4x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t13 = sin(pkin(4));
t33 = g(3) * t13;
t17 = sin(qJ(3));
t32 = t13 * t17;
t18 = sin(qJ(2));
t31 = t13 * t18;
t20 = cos(qJ(3));
t30 = t13 * t20;
t21 = cos(qJ(2));
t29 = t13 * t21;
t15 = cos(pkin(4));
t28 = t15 * t18;
t27 = t15 * t21;
t16 = sin(qJ(4));
t26 = t16 * t20;
t19 = cos(qJ(4));
t25 = t19 * t20;
t24 = t20 * t21;
t12 = sin(pkin(8));
t14 = cos(pkin(8));
t6 = t12 * t21 + t14 * t28;
t8 = -t12 * t28 + t14 * t21;
t23 = g(1) * (t12 * t30 - t8 * t17) + g(2) * (-t14 * t30 - t6 * t17) + g(3) * (t15 * t20 - t17 * t31);
t5 = t12 * t18 - t14 * t27;
t7 = t12 * t27 + t14 * t18;
t22 = -g(1) * t7 - g(2) * t5 + g(3) * t29;
t10 = t15 * t17 + t18 * t30;
t4 = t12 * t32 + t8 * t20;
t2 = -t14 * t32 + t6 * t20;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t22, g(1) * t8 + g(2) * t6 + g(3) * t31, 0, 0, 0, 0, 0, -t22 * t20, t22 * t17, 0, 0, 0, 0, 0, -g(1) * (t8 * t16 - t7 * t25) - g(2) * (t6 * t16 - t5 * t25) - (t16 * t18 + t19 * t24) * t33, -g(1) * (t8 * t19 + t7 * t26) - g(2) * (t6 * t19 + t5 * t26) - (-t16 * t24 + t18 * t19) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, g(1) * t4 + g(2) * t2 + g(3) * t10, 0, 0, 0, 0, 0, -t23 * t19, t23 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t4 * t16 + t7 * t19) - g(2) * (-t2 * t16 + t5 * t19) - g(3) * (-t10 * t16 - t19 * t29), -g(1) * (-t7 * t16 - t4 * t19) - g(2) * (-t5 * t16 - t2 * t19) - g(3) * (-t10 * t19 + t16 * t29);];
taug_reg = t1;
