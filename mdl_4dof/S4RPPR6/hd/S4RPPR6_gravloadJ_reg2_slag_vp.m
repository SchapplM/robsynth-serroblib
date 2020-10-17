% Calculate inertial parameters regressor of gravitation load for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:44
% EndTime: 2019-12-31 16:40:44
% DurationCPUTime: 0.12s
% Computational Cost: add. (53->32), mult. (128->48), div. (0->0), fcn. (132->6), ass. (0->25)
t19 = sin(qJ(1));
t29 = g(1) * t19;
t17 = cos(pkin(6));
t21 = cos(qJ(1));
t28 = t17 * t21;
t27 = t21 * pkin(1) + t19 * qJ(2);
t16 = sin(pkin(6));
t26 = qJ(3) * t16;
t25 = pkin(2) * t28 + t21 * t26 + t27;
t9 = g(1) * t21 + g(2) * t19;
t8 = -g(2) * t21 + t29;
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t24 = t16 * t20 - t17 * t18;
t23 = t16 * t18 + t17 * t20;
t22 = -pkin(2) * t17 - pkin(1) - t26;
t13 = t21 * qJ(2);
t7 = t8 * t17;
t6 = t8 * t16;
t5 = t23 * t21;
t4 = t24 * t21;
t3 = t23 * t19;
t2 = t24 * t19;
t1 = g(3) * t17 - t9 * t16;
t10 = [0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, -t9, -g(1) * (-t19 * pkin(1) + t13) - g(2) * t27, 0, 0, 0, 0, 0, 0, t7, -t9, t6, -g(1) * t13 - g(2) * t25 - t22 * t29, 0, 0, 0, 0, 0, 0, g(1) * t3 - g(2) * t5, g(1) * t2 - g(2) * t4, t9, -g(1) * (-t21 * pkin(5) + t13) - g(2) * (pkin(3) * t28 + t25) + (-g(1) * (-pkin(3) * t17 + t22) + g(2) * pkin(5)) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2 + g(3) * t23, g(1) * t5 + g(2) * t3 + g(3) * t24, 0, 0;];
taug_reg = t10;
