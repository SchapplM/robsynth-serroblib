% Calculate inertial parameters regressor of gravitation load for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:55
% EndTime: 2019-12-05 15:26:56
% DurationCPUTime: 0.15s
% Computational Cost: add. (107->39), mult. (136->48), div. (0->0), fcn. (131->8), ass. (0->26)
t13 = sin(pkin(7));
t14 = cos(pkin(7));
t22 = g(1) * t14 + g(2) * t13;
t12 = qJ(2) + pkin(8);
t10 = cos(t12);
t9 = sin(t12);
t2 = g(3) * t9 + t22 * t10;
t16 = sin(qJ(2));
t33 = pkin(2) * t16;
t30 = g(3) * t10;
t15 = sin(qJ(5));
t29 = t13 * t15;
t17 = cos(qJ(5));
t28 = t13 * t17;
t27 = t14 * t15;
t26 = t14 * t17;
t25 = qJ(4) * t10;
t18 = cos(qJ(2));
t24 = t18 * pkin(2) + t10 * pkin(3) + t9 * qJ(4);
t23 = -pkin(3) * t9 - t33;
t19 = -g(3) * t18 + t22 * t16;
t5 = t14 * t25;
t4 = t13 * t25;
t3 = -g(1) * t13 + g(2) * t14;
t1 = t22 * t9 - t30;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, g(3) * t16 + t22 * t18, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t19 * pkin(2), 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (t23 * t14 + t5) - g(2) * (t23 * t13 + t4) - g(3) * t24, 0, 0, 0, 0, 0, 0, -t2 * t15, -t2 * t17, t1, -g(1) * t5 - g(2) * t4 - g(3) * (t10 * pkin(6) + t24) + t22 * (t33 - (-pkin(3) - pkin(6)) * t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t9 * t26 - t29) - g(2) * (t9 * t28 + t27) + t17 * t30, -g(1) * (-t9 * t27 - t28) - g(2) * (-t9 * t29 + t26) - t15 * t30, 0, 0;];
taug_reg = t6;
