% Calculate inertial parameters regressor of gravitation load for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:12
% EndTime: 2019-12-31 18:11:12
% DurationCPUTime: 0.16s
% Computational Cost: add. (149->44), mult. (161->52), div. (0->0), fcn. (147->6), ass. (0->29)
t21 = sin(qJ(3));
t20 = qJ(1) + pkin(7);
t14 = sin(t20);
t15 = cos(t20);
t6 = g(1) * t15 + g(2) * t14;
t39 = t6 * t21;
t16 = t21 * qJ(4);
t23 = cos(qJ(3));
t31 = t23 * pkin(3) + t16;
t37 = pkin(3) * t21;
t36 = g(1) * t14;
t33 = t23 * pkin(4);
t32 = t15 * t23;
t30 = qJ(4) * t23;
t24 = cos(qJ(1));
t29 = t24 * pkin(1) + t15 * pkin(2) + t14 * pkin(6);
t22 = sin(qJ(1));
t28 = -t22 * pkin(1) + t15 * pkin(6);
t27 = pkin(3) * t32 + t15 * t16 + t29;
t5 = -g(2) * t15 + t36;
t26 = g(1) * t22 - g(2) * t24;
t25 = -pkin(2) - t31;
t9 = t15 * t30;
t7 = t14 * t30;
t4 = t5 * t23;
t3 = t5 * t21;
t2 = g(3) * t21 + t6 * t23;
t1 = -g(3) * t23 + t39;
t8 = [0, 0, 0, 0, 0, 0, t26, g(1) * t24 + g(2) * t22, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t26 * pkin(1), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t14 * pkin(2) + t28) - g(2) * t29, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t28 - g(2) * t27 - t25 * t36, 0, 0, 0, 0, 0, 0, t4, t3, t6, -g(1) * (-t15 * qJ(5) + t28) - g(2) * (pkin(4) * t32 + t27) + (-g(1) * (t25 - t33) + g(2) * qJ(5)) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-t15 * t37 + t9) - g(2) * (-t14 * t37 + t7) - g(3) * t31, 0, 0, 0, 0, 0, 0, t1, -t2, 0, -g(1) * t9 - g(2) * t7 - g(3) * (t31 + t33) + (pkin(3) + pkin(4)) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5;];
taug_reg = t8;
