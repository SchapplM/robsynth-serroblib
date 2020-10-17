% Calculate inertial parameters regressor of gravitation load for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:41
% EndTime: 2019-12-05 15:33:42
% DurationCPUTime: 0.15s
% Computational Cost: add. (128->37), mult. (175->45), div. (0->0), fcn. (171->8), ass. (0->28)
t13 = sin(pkin(7));
t14 = cos(pkin(7));
t23 = g(1) * t14 + g(2) * t13;
t12 = qJ(2) + pkin(8);
t10 = cos(t12);
t9 = sin(t12);
t5 = -g(3) * t10 + t23 * t9;
t32 = g(3) * t9;
t17 = sin(qJ(2));
t31 = pkin(2) * t17;
t16 = sin(qJ(4));
t27 = t13 * t16;
t18 = cos(qJ(4));
t26 = t13 * t18;
t25 = t14 * t16;
t24 = t14 * t18;
t19 = cos(qJ(2));
t20 = -g(3) * t19 + t23 * t17;
t1 = -g(1) * (-t10 * t25 + t26) - g(2) * (-t10 * t27 - t24) + t16 * t32;
t15 = -qJ(5) - pkin(6);
t11 = t19 * pkin(2);
t8 = t18 * pkin(4) + pkin(3);
t7 = -g(1) * t13 + g(2) * t14;
t6 = t23 * t10 + t32;
t4 = t5 * t18;
t3 = t5 * t16;
t2 = -g(1) * (-t10 * t24 - t27) - g(2) * (-t10 * t26 + t25) + t18 * t32;
t21 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, g(3) * t17 + t23 * t19, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t20 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * (t10 * pkin(3) + t9 * pkin(6) + t11) + t23 * (pkin(3) * t9 - pkin(6) * t10 + t31), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * (t10 * t8 - t9 * t15 + t11) + t23 * (t10 * t15 + t8 * t9 + t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t21;
