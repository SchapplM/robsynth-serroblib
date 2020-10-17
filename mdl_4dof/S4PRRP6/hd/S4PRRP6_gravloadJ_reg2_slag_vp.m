% Calculate inertial parameters regressor of gravitation load for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:43
% EndTime: 2019-12-31 16:30:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (65->36), mult. (167->51), div. (0->0), fcn. (172->6), ass. (0->28)
t18 = sin(qJ(2));
t15 = sin(pkin(6));
t16 = cos(pkin(6));
t25 = g(1) * t16 + g(2) * t15;
t34 = t25 * t18;
t33 = pkin(2) * t18;
t20 = cos(qJ(2));
t32 = pkin(5) * t20;
t29 = g(3) * t18;
t17 = sin(qJ(3));
t28 = t17 * t20;
t19 = cos(qJ(3));
t27 = t19 * t20;
t26 = t20 * pkin(2) + t18 * pkin(5);
t24 = pkin(3) * t19 + qJ(4) * t17;
t5 = t15 * t28 + t16 * t19;
t7 = -t15 * t19 + t16 * t28;
t1 = g(1) * t7 + g(2) * t5 + t17 * t29;
t6 = t15 * t27 - t16 * t17;
t8 = t15 * t17 + t16 * t27;
t22 = g(1) * t8 + g(2) * t6 + t19 * t29;
t21 = -g(3) * t20 + t34;
t10 = t16 * t32;
t9 = t15 * t32;
t4 = t25 * t20 + t29;
t3 = t21 * t19;
t2 = t21 * t17;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t2, -t4, -g(1) * (-t16 * t33 + t10) - g(2) * (-t15 * t33 + t9) - g(3) * t26, 0, 0, 0, 0, 0, 0, t3, -t4, t2, -g(1) * t10 - g(2) * t9 - g(3) * (t24 * t20 + t26) + (pkin(2) + t24) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t22, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t22, -g(1) * (-t7 * pkin(3) + t8 * qJ(4)) - g(2) * (-t5 * pkin(3) + t6 * qJ(4)) - (-pkin(3) * t17 + qJ(4) * t19) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t11;
