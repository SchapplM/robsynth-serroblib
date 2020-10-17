% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR15_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:08
% EndTime: 2019-12-31 20:43:09
% DurationCPUTime: 0.18s
% Computational Cost: add. (107->42), mult. (198->62), div. (0->0), fcn. (213->8), ass. (0->37)
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t15 = g(1) * t25 + g(2) * t22;
t21 = sin(qJ(2));
t24 = cos(qJ(2));
t8 = g(3) * t21 + t15 * t24;
t39 = g(3) * t24;
t19 = qJ(4) + qJ(5);
t16 = sin(t19);
t38 = t22 * t16;
t17 = cos(t19);
t37 = t22 * t17;
t20 = sin(qJ(4));
t36 = t22 * t20;
t23 = cos(qJ(4));
t35 = t22 * t23;
t34 = t25 * t16;
t33 = t25 * t17;
t32 = t25 * t20;
t31 = t25 * t23;
t30 = g(1) * t22 - g(2) * t25;
t29 = t24 * pkin(2) + t21 * qJ(3);
t27 = pkin(1) + t29;
t14 = t30 * t24;
t13 = t30 * t21;
t12 = -t21 * t36 + t31;
t11 = t21 * t35 + t32;
t10 = t21 * t32 + t35;
t9 = t21 * t31 - t36;
t7 = t15 * t21 - t39;
t6 = -t21 * t38 + t33;
t5 = t21 * t37 + t34;
t4 = t21 * t34 + t37;
t3 = t21 * t33 - t38;
t2 = g(1) * t4 - g(2) * t6 - t16 * t39;
t1 = -g(1) * t3 - g(2) * t5 + t17 * t39;
t18 = [0, t30, t15, 0, 0, 0, 0, 0, t14, -t13, -t15, -t14, t13, (-g(1) * pkin(6) - g(2) * t27) * t25 + (-g(2) * pkin(6) + g(1) * t27) * t22, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, -t7, -t8, -g(3) * t29 + t15 * (pkin(2) * t21 - qJ(3) * t24), 0, 0, 0, 0, 0, -t8 * t20, -t8 * t23, 0, 0, 0, 0, 0, -t8 * t16, -t8 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11 + t23 * t39, g(1) * t10 - g(2) * t12 - t20 * t39, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t18;
