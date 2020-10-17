% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:21
% EndTime: 2019-12-31 18:07:22
% DurationCPUTime: 0.18s
% Computational Cost: add. (118->48), mult. (172->59), div. (0->0), fcn. (169->8), ass. (0->31)
t23 = sin(qJ(1));
t25 = cos(qJ(1));
t43 = -g(1) * t23 + g(2) * t25;
t18 = pkin(8) + qJ(4);
t12 = sin(t18);
t13 = cos(t18);
t42 = -g(3) * t12 - t13 * t43;
t19 = sin(pkin(8));
t41 = pkin(3) * t19;
t37 = g(3) * t13;
t22 = sin(qJ(5));
t36 = t23 * t22;
t24 = cos(qJ(5));
t35 = t23 * t24;
t34 = t25 * t22;
t33 = t25 * t24;
t32 = t25 * pkin(1) + t23 * qJ(2);
t15 = t25 * qJ(2);
t31 = -t23 * pkin(1) + t15;
t29 = t12 * pkin(4) - t13 * pkin(7);
t8 = g(1) * t25 + g(2) * t23;
t21 = -pkin(6) - qJ(3);
t28 = t23 * t21 + t25 * t41 + t31;
t27 = -t25 * t21 + t23 * t41 + t32;
t6 = t12 * t33 - t36;
t5 = t12 * t34 + t35;
t4 = t12 * t35 + t34;
t3 = -t12 * t36 + t33;
t2 = t8 * t13;
t1 = -t12 * t43 + t37;
t7 = [0, 0, 0, 0, 0, 0, -t43, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t8, -g(1) * t31 - g(2) * t32, 0, 0, 0, 0, 0, 0, -t8 * t19, -t8 * cos(pkin(8)), -t43, -g(1) * (t15 + (-pkin(1) - qJ(3)) * t23) - g(2) * (t25 * qJ(3) + t32), 0, 0, 0, 0, 0, 0, -t8 * t12, -t2, -t43, -g(1) * t28 - g(2) * t27, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t2, -g(1) * (t29 * t25 + t28) - g(2) * (t29 * t23 + t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t42 * t24, t42 * t22, -t1, g(3) * t29 + t43 * (pkin(4) * t13 + pkin(7) * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t22 * t37, g(1) * t4 - g(2) * t6 + t24 * t37, 0, 0;];
taug_reg = t7;
