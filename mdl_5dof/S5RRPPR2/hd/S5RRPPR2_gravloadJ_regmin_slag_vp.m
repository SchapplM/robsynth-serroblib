% Calculate minimal parameter regressor of gravitation load for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:36
% EndTime: 2020-01-03 11:57:36
% DurationCPUTime: 0.13s
% Computational Cost: add. (159->33), mult. (114->50), div. (0->0), fcn. (118->10), ass. (0->34)
t25 = sin(pkin(9));
t36 = g(1) * t25;
t26 = cos(pkin(9));
t27 = sin(qJ(5));
t35 = t26 * t27;
t29 = cos(qJ(5));
t34 = t26 * t29;
t24 = qJ(1) + qJ(2);
t19 = pkin(8) + t24;
t15 = sin(t19);
t16 = cos(t19);
t21 = cos(t24);
t18 = pkin(2) * t21;
t33 = t16 * pkin(3) + t15 * qJ(4) + t18;
t20 = sin(t24);
t17 = pkin(2) * t20;
t32 = t15 * pkin(3) - t16 * qJ(4) + t17;
t31 = g(2) * t16 + g(3) * t15;
t11 = -g(2) * t21 - g(3) * t20;
t30 = cos(qJ(1));
t28 = sin(qJ(1));
t23 = t30 * pkin(1);
t22 = t28 * pkin(1);
t10 = g(2) * t20 - g(3) * t21;
t9 = -g(2) * t15 + g(3) * t16;
t8 = t31 * t26;
t7 = t31 * t25;
t6 = t15 * t27 + t16 * t34;
t5 = -t15 * t29 + t16 * t35;
t4 = t15 * t34 - t16 * t27;
t3 = -t15 * t35 - t16 * t29;
t2 = -g(2) * t6 - g(3) * t4;
t1 = g(2) * t5 - g(3) * t3;
t12 = [0, -g(2) * t30 - g(3) * t28, g(2) * t28 - g(3) * t30, 0, t11, t10, -g(2) * (t18 + t23) - g(3) * (t17 + t22), -t8, t7, t9, -g(2) * (t23 + t33) - g(3) * (t22 + t32), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, t11, t10, t11 * pkin(2), -t8, t7, t9, -g(2) * t33 - g(3) * t32, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t3 - g(3) * t5 + t27 * t36, g(2) * t4 - g(3) * t6 + t29 * t36;];
taug_reg = t12;
