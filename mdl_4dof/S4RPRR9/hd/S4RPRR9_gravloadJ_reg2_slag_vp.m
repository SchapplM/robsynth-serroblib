% Calculate inertial parameters regressor of gravitation load for
% S4RPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR9_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR9_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:56:29
% EndTime: 2019-12-31 16:56:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (57->36), mult. (138->51), div. (0->0), fcn. (137->6), ass. (0->26)
t15 = sin(qJ(3));
t18 = cos(qJ(3));
t21 = t15 * pkin(3) - t18 * pkin(6);
t16 = sin(qJ(1));
t19 = cos(qJ(1));
t37 = -g(1) * t16 + g(2) * t19;
t36 = -g(3) * t15 - t18 * t37;
t35 = -pkin(1) - pkin(5);
t31 = g(3) * t18;
t28 = t19 * pkin(1) + t16 * qJ(2);
t14 = sin(qJ(4));
t27 = t16 * t14;
t17 = cos(qJ(4));
t26 = t16 * t17;
t25 = t19 * t14;
t24 = t19 * t17;
t23 = g(2) * (t19 * pkin(5) + t28);
t8 = g(1) * t19 + g(2) * t16;
t10 = t19 * qJ(2);
t6 = t8 * t18;
t5 = t15 * t24 - t27;
t4 = t15 * t25 + t26;
t3 = t15 * t26 + t25;
t2 = -t15 * t27 + t24;
t1 = -t15 * t37 + t31;
t7 = [0, 0, 0, 0, 0, 0, -t37, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t8, -g(1) * (-t16 * pkin(1) + t10) - g(2) * t28, 0, 0, 0, 0, 0, 0, -t8 * t15, -t6, -t37, -g(1) * (t35 * t16 + t10) - t23, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t3, g(1) * t4 - g(2) * t2, t6, -g(1) * (t21 * t19 + t10) - t23 + (-g(1) * t35 - g(2) * t21) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t36 * t17, t36 * t14, -t1, g(3) * t21 + t37 * (pkin(3) * t18 + pkin(6) * t15); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4 + t14 * t31, g(1) * t3 - g(2) * t5 + t17 * t31, 0, 0;];
taug_reg = t7;
