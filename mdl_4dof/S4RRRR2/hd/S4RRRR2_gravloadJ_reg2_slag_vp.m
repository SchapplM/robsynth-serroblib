% Calculate inertial parameters regressor of gravitation load for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:18
% EndTime: 2019-12-31 17:23:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (139->32), mult. (126->41), div. (0->0), fcn. (113->8), ass. (0->29)
t22 = sin(qJ(1));
t32 = t22 * pkin(1);
t20 = qJ(1) + qJ(2);
t15 = sin(t20);
t17 = cos(t20);
t31 = t17 * pkin(2) + t15 * pkin(6);
t30 = -t15 * pkin(2) + t17 * pkin(6);
t23 = cos(qJ(3));
t13 = t23 * pkin(3) + pkin(2);
t25 = -pkin(7) - pkin(6);
t29 = t17 * t13 - t15 * t25;
t8 = g(1) * t17 + g(2) * t15;
t7 = g(1) * t15 - g(2) * t17;
t24 = cos(qJ(1));
t28 = g(1) * t22 - g(2) * t24;
t27 = -t15 * t13 - t17 * t25;
t21 = sin(qJ(3));
t26 = -g(3) * t23 + t8 * t21;
t19 = qJ(3) + qJ(4);
t18 = t24 * pkin(1);
t16 = cos(t19);
t14 = sin(t19);
t6 = t7 * t23;
t5 = t7 * t21;
t4 = t7 * t16;
t3 = t7 * t14;
t2 = g(3) * t14 + t8 * t16;
t1 = -g(3) * t16 + t8 * t14;
t9 = [0, 0, 0, 0, 0, 0, t28, g(1) * t24 + g(2) * t22, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t28 * pkin(1), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t30 - t32) - g(2) * (t18 + t31), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t27 - t32) - g(2) * (t18 + t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t30 - g(2) * t31, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * t27 - g(2) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, g(3) * t21 + t8 * t23, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t26 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t9;
