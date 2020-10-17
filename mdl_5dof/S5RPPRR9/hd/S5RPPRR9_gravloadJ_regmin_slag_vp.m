% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR9
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
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:48
% EndTime: 2019-12-31 18:02:48
% DurationCPUTime: 0.12s
% Computational Cost: add. (72->31), mult. (164->50), div. (0->0), fcn. (200->8), ass. (0->25)
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t29 = sin(qJ(1));
t30 = cos(qJ(1));
t3 = -t29 * t24 - t30 * t25;
t4 = t30 * t24 - t29 * t25;
t22 = g(1) * t3 + g(2) * t4;
t33 = -g(3) * t17 + t22 * t15;
t32 = g(3) * t15;
t14 = sin(qJ(5));
t28 = t14 * t17;
t16 = cos(qJ(5));
t27 = t16 * t17;
t26 = t30 * pkin(1) + t29 * qJ(2);
t23 = g(1) * t4 - g(2) * t3;
t21 = -t29 * pkin(1) + t30 * qJ(2);
t20 = t3 * t14 + t4 * t27;
t19 = -t3 * t16 + t4 * t28;
t6 = g(1) * t30 + g(2) * t29;
t5 = g(1) * t29 - g(2) * t30;
t2 = t4 * t14 - t3 * t27;
t1 = t4 * t16 + t3 * t28;
t7 = [0, t5, t6, t5, -t6, -g(1) * t21 - g(2) * t26, -t23, t22, -g(1) * (-t29 * pkin(2) + t21) - g(2) * (t30 * pkin(2) + t26), 0, 0, 0, 0, 0, -t23 * t17, t23 * t15, 0, 0, 0, 0, 0, -g(1) * t20 - g(2) * t2, g(1) * t19 - g(2) * t1; 0, 0, 0, 0, 0, -t5, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t22 * t17 - t32, 0, 0, 0, 0, 0, -t33 * t16, t33 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t19 - t14 * t32, g(1) * t2 - g(2) * t20 - t16 * t32;];
taug_reg = t7;
