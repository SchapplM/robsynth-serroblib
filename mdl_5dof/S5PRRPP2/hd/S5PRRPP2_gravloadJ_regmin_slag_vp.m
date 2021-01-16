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
% taug_reg [5x19]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% OptimizationMode: 2
% StartTime: 2021-01-15 15:32:32
% EndTime: 2021-01-15 15:32:34
% DurationCPUTime: 0.20s
% Computational Cost: add. (154->56), mult. (246->80), div. (0->0), fcn. (252->8), ass. (0->39)
t20 = qJ(3) + pkin(8);
t17 = sin(t20);
t18 = cos(t20);
t47 = pkin(4) * t18 + qJ(5) * t17;
t21 = sin(pkin(7));
t22 = cos(pkin(7));
t31 = g(1) * t22 + g(2) * t21;
t25 = sin(qJ(2));
t27 = cos(qJ(2));
t8 = -g(3) * t27 + t31 * t25;
t43 = g(3) * t25;
t26 = cos(qJ(3));
t16 = t26 * pkin(3) + pkin(2);
t41 = t16 * t25;
t40 = t21 * t26;
t39 = t21 * t27;
t38 = t22 * t26;
t37 = t22 * t27;
t23 = qJ(4) + pkin(6);
t36 = t23 * t27;
t24 = sin(qJ(3));
t35 = t24 * t27;
t34 = t26 * t27;
t32 = t22 * t35;
t30 = -t21 * t35 - t38;
t4 = t17 * t39 + t22 * t18;
t6 = t17 * t37 - t21 * t18;
t1 = g(1) * t6 + g(2) * t4 + t17 * t43;
t5 = -t22 * t17 + t18 * t39;
t7 = t21 * t17 + t18 * t37;
t28 = g(1) * t7 + g(2) * t5 + t18 * t43;
t9 = t31 * t27 + t43;
t15 = pkin(3) * t40;
t12 = t27 * t16;
t11 = t22 * t36;
t10 = t21 * t36;
t3 = t8 * t18;
t2 = t8 * t17;
t13 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, t8, t9, 0, 0, 0, 0, 0, t8 * t26, -t8 * t24, t3, -t2, -t9, -g(1) * (-t22 * t41 + t11) - g(2) * (-t21 * t41 + t10) - g(3) * (t25 * t23 + t12), t3, -t9, t2, -g(1) * t11 - g(2) * t10 - g(3) * (t47 * t27 + t12) + (-g(3) * t23 + t31 * (t16 + t47)) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t32 + t40) - g(2) * t30 + t24 * t43, -g(1) * (-t21 * t24 - t22 * t34) - g(2) * (-t21 * t34 + t22 * t24) + t26 * t43, t1, t28, 0, -g(1) * t15 + (g(2) * t38 + t9 * t24) * pkin(3), t1, 0, -t28, -g(1) * (-pkin(3) * t32 - t6 * pkin(4) + t7 * qJ(5) + t15) - g(2) * (t30 * pkin(3) - t4 * pkin(4) + t5 * qJ(5)) - (-pkin(3) * t24 - pkin(4) * t17 + qJ(5) * t18) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t13;
