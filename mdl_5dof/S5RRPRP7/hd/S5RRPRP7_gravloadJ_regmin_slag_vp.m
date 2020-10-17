% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:31
% EndTime: 2019-12-31 20:01:32
% DurationCPUTime: 0.19s
% Computational Cost: add. (171->52), mult. (256->72), div. (0->0), fcn. (269->8), ass. (0->37)
t19 = qJ(2) + pkin(8);
t16 = cos(t19);
t15 = sin(t19);
t39 = t15 * pkin(7);
t43 = t16 * pkin(3) + t39;
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t10 = g(1) * t26 + g(2) * t23;
t40 = g(3) * t15;
t21 = sin(qJ(4));
t37 = t23 * t21;
t24 = cos(qJ(4));
t36 = t23 * t24;
t20 = -qJ(3) - pkin(6);
t35 = t26 * t20;
t34 = t26 * t21;
t33 = t26 * t24;
t5 = t16 * t37 + t33;
t7 = t16 * t34 - t36;
t32 = g(1) * t5 - g(2) * t7;
t9 = g(1) * t23 - g(2) * t26;
t31 = pkin(4) * t24 + qJ(5) * t21 + pkin(3);
t1 = g(1) * t7 + g(2) * t5 + t21 * t40;
t6 = t16 * t36 - t34;
t8 = t16 * t33 + t37;
t30 = g(1) * t8 + g(2) * t6 + t24 * t40;
t29 = -g(3) * t16 + t10 * t15;
t22 = sin(qJ(2));
t25 = cos(qJ(2));
t28 = -g(3) * t25 + t10 * t22;
t17 = t25 * pkin(2);
t14 = t17 + pkin(1);
t11 = t26 * t14;
t4 = t29 * t24;
t3 = t29 * t21;
t2 = g(1) * t6 - g(2) * t8;
t12 = [0, t9, t10, 0, 0, 0, 0, 0, t9 * t25, -t9 * t22, -t10, -g(1) * (-t23 * t14 - t35) - g(2) * (-t23 * t20 + t11), 0, 0, 0, 0, 0, t2, -t32, t2, t9 * t15, t32, -g(1) * (-t6 * pkin(4) - t5 * qJ(5) - t35) - g(2) * (t8 * pkin(4) + t7 * qJ(5) + t43 * t26 + t11) + (-g(1) * (-t14 - t43) + g(2) * t20) * t23; 0, 0, 0, 0, 0, 0, 0, 0, t28, g(3) * t22 + t10 * t25, 0, t28 * pkin(2), 0, 0, 0, 0, 0, t4, -t3, t4, -t10 * t16 - t40, t3, -g(3) * (t16 * t31 + t17 + t39) + t10 * (pkin(2) * t22 - pkin(7) * t16 + t15 * t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t30, t1, 0, -t30, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (-t5 * pkin(4) + t6 * qJ(5)) - (-pkin(4) * t21 + qJ(5) * t24) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t12;
