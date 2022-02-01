% Calculate inertial parameters regressor of gravitation load for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:45
% EndTime: 2022-01-23 09:12:46
% DurationCPUTime: 0.20s
% Computational Cost: add. (160->44), mult. (174->60), div. (0->0), fcn. (179->8), ass. (0->31)
t21 = sin(pkin(8));
t22 = cos(pkin(8));
t26 = cos(qJ(4));
t42 = -t21 * (-qJ(5) - pkin(6)) + (t26 * pkin(4) + pkin(3)) * t22;
t24 = sin(qJ(4));
t38 = g(3) * t21;
t20 = qJ(1) + pkin(7);
t17 = sin(t20);
t18 = cos(t20);
t34 = t22 * t24;
t5 = t17 * t34 + t18 * t26;
t7 = t17 * t26 - t18 * t34;
t1 = -g(1) * t7 + g(2) * t5 + t24 * t38;
t39 = g(1) * t17;
t36 = t18 * t24;
t33 = t22 * t26;
t27 = cos(qJ(1));
t31 = t27 * pkin(1) + t18 * pkin(2) + t17 * qJ(3);
t25 = sin(qJ(1));
t30 = -t25 * pkin(1) + t18 * qJ(3);
t29 = pkin(3) * t22 + pkin(6) * t21;
t11 = g(1) * t18 + g(2) * t17;
t10 = -g(2) * t18 + t39;
t28 = g(1) * t25 - g(2) * t27;
t9 = t10 * t21;
t8 = t17 * t24 + t18 * t33;
t6 = -t17 * t33 + t36;
t4 = -g(1) * t6 - g(2) * t8;
t3 = -g(1) * t5 - g(2) * t7;
t2 = g(1) * t8 - g(2) * t6 + t26 * t38;
t12 = [0, 0, 0, 0, 0, 0, t28, g(1) * t27 + g(2) * t25, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, t28 * pkin(1), 0, 0, 0, 0, 0, 0, t10 * t22, -t9, -t11, -g(1) * (-t17 * pkin(2) + t30) - g(2) * t31, 0, 0, 0, 0, 0, 0, t4, t3, t9, -g(1) * t30 - g(2) * (t18 * t29 + t31) - (-pkin(2) - t29) * t39, 0, 0, 0, 0, 0, 0, t4, t3, t9, -g(1) * (pkin(4) * t36 + t30) - g(2) * (t42 * t18 + t31) + (-g(1) * (-pkin(2) - t42) - g(2) * pkin(4) * t24) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t22 - t11 * t21;];
taug_reg = t12;
