% Calculate minimal parameter regressor of gravitation load for
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
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:55:47
% EndTime: 2021-01-15 16:55:48
% DurationCPUTime: 0.13s
% Computational Cost: add. (127->35), mult. (144->54), div. (0->0), fcn. (151->8), ass. (0->28)
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t24 = cos(qJ(4));
t39 = -t19 * (-qJ(5) - pkin(6)) + (t24 * pkin(4) + pkin(3)) * t20;
t22 = sin(qJ(4));
t36 = g(1) * t19;
t18 = qJ(1) + pkin(7);
t14 = sin(t18);
t15 = cos(t18);
t32 = t20 * t22;
t5 = -t14 * t32 - t15 * t24;
t7 = -t14 * t24 + t15 * t32;
t1 = -g(2) * t5 - g(3) * t7 + t22 * t36;
t34 = t14 * t22;
t31 = t20 * t24;
t23 = sin(qJ(1));
t30 = t23 * pkin(1) + t14 * pkin(2);
t25 = cos(qJ(1));
t28 = t25 * pkin(1) + t15 * pkin(2) + t14 * qJ(3);
t9 = g(2) * t15 + g(3) * t14;
t27 = -g(2) * t14 + g(3) * t15;
t26 = -g(2) * t25 - g(3) * t23;
t8 = t15 * t31 + t34;
t6 = t14 * t31 - t15 * t22;
t4 = -g(2) * t8 - g(3) * t6;
t3 = g(2) * t7 - g(3) * t5;
t2 = g(2) * t6 - g(3) * t8 + t24 * t36;
t10 = [0, t26, g(2) * t23 - g(3) * t25, t26 * pkin(1), -t9 * t20, t27, -g(2) * t28 - g(3) * (-t15 * qJ(3) + t30), 0, 0, 0, 0, 0, t4, t3, t4, t3, -t9 * t19, -g(2) * (pkin(4) * t34 + t28) - g(3) * (t39 * t14 + t30) + (-g(2) * t39 - g(3) * (-pkin(4) * t22 - qJ(3))) * t15; 0, 0, 0, -g(1), 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t20 + t27 * t19;];
taug_reg = t10;
