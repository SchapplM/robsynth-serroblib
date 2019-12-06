% Calculate minimal parameter regressor of gravitation load for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = sin(pkin(7));
t16 = cos(pkin(7));
t27 = g(1) * t16 + g(2) * t15;
t19 = sin(qJ(2));
t21 = cos(qJ(2));
t5 = -g(3) * t21 + t27 * t19;
t38 = g(3) * t19;
t20 = cos(qJ(3));
t36 = t15 * t20;
t35 = t15 * t21;
t34 = t16 * t20;
t33 = t16 * t21;
t17 = -qJ(4) - pkin(6);
t32 = t17 * t21;
t18 = sin(qJ(3));
t31 = t18 * t21;
t30 = t20 * t21;
t29 = t16 * t31;
t10 = t20 * pkin(3) + pkin(2);
t28 = t21 * t10 - t19 * t17;
t14 = qJ(3) + pkin(8);
t11 = sin(t14);
t12 = cos(t14);
t26 = pkin(4) * t12 + qJ(5) * t11;
t24 = -t15 * t31 - t34;
t1 = t11 * t35 + t16 * t12;
t3 = t11 * t33 - t15 * t12;
t23 = g(1) * t3 + g(2) * t1 + t11 * t38;
t6 = t27 * t21 + t38;
t9 = pkin(3) * t36;
t4 = t15 * t11 + t12 * t33;
t2 = -t16 * t11 + t12 * t35;
t7 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, t5, t6, 0, 0, 0, 0, 0, t5 * t20, -t5 * t18, -t6, -g(3) * t28 + t27 * (t10 * t19 + t32), t5 * t12, -t6, t5 * t11, -g(3) * (t26 * t21 + t28) + t27 * (t32 - (-t10 - t26) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t29 + t36) - g(2) * t24 + t18 * t38, -g(1) * (-t15 * t18 - t16 * t30) - g(2) * (-t15 * t30 + t16 * t18) + t20 * t38, 0, -g(1) * t9 + (g(2) * t34 + t6 * t18) * pkin(3), t23, 0, -g(1) * t4 - g(2) * t2 - t12 * t38, -g(1) * (-pkin(3) * t29 - t3 * pkin(4) + t4 * qJ(5) + t9) - g(2) * (t24 * pkin(3) - t1 * pkin(4) + t2 * qJ(5)) - (-pkin(3) * t18 - pkin(4) * t11 + qJ(5) * t12) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23;];
taug_reg = t7;
