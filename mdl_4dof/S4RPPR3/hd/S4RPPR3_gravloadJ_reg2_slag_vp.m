% Calculate inertial parameters regressor of gravitation load for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:56
% EndTime: 2019-12-31 16:37:57
% DurationCPUTime: 0.09s
% Computational Cost: add. (76->29), mult. (64->32), div. (0->0), fcn. (58->8), ass. (0->17)
t15 = sin(qJ(1));
t18 = t15 * pkin(1);
t11 = qJ(1) + pkin(6);
t6 = sin(t11);
t8 = cos(t11);
t2 = g(1) * t8 + g(2) * t6;
t1 = g(1) * t6 - g(2) * t8;
t16 = cos(qJ(1));
t17 = g(1) * t15 - g(2) * t16;
t14 = -pkin(5) - qJ(3);
t13 = cos(pkin(7));
t10 = pkin(7) + qJ(4);
t9 = t16 * pkin(1);
t7 = cos(t10);
t5 = sin(t10);
t4 = t13 * pkin(3) + pkin(2);
t3 = [0, 0, 0, 0, 0, 0, t17, g(1) * t16 + g(2) * t15, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t17 * pkin(1), 0, 0, 0, 0, 0, 0, t1 * t13, -t1 * sin(pkin(7)), -t2, -g(1) * (-t6 * pkin(2) + t8 * qJ(3) - t18) - g(2) * (t8 * pkin(2) + t6 * qJ(3) + t9), 0, 0, 0, 0, 0, 0, t1 * t7, -t1 * t5, -t2, -g(1) * (-t8 * t14 - t6 * t4 - t18) - g(2) * (-t6 * t14 + t8 * t4 + t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t7 + t2 * t5, g(3) * t5 + t2 * t7, 0, 0;];
taug_reg = t3;
