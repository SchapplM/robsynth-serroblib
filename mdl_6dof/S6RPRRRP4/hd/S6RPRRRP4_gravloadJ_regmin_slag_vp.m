% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:25:29
% EndTime: 2019-05-06 01:25:30
% DurationCPUTime: 0.20s
% Computational Cost: add. (256->48), mult. (233->64), div. (0->0), fcn. (228->10), ass. (0->36)
t23 = pkin(10) + qJ(3);
t19 = cos(t23);
t20 = qJ(4) + t23;
t15 = sin(t20);
t16 = cos(t20);
t29 = cos(qJ(5));
t17 = t29 * pkin(5) + pkin(4);
t26 = -qJ(6) - pkin(9);
t34 = -t15 * t26 + t16 * t17;
t49 = pkin(3) * t19 + t34;
t28 = sin(qJ(1));
t30 = cos(qJ(1));
t13 = g(1) * t30 + g(2) * t28;
t27 = sin(qJ(5));
t43 = g(3) * t15;
t37 = t30 * t29;
t40 = t28 * t27;
t6 = t16 * t40 + t37;
t38 = t30 * t27;
t39 = t28 * t29;
t8 = -t16 * t38 + t39;
t48 = -g(1) * t8 + g(2) * t6 + t27 * t43;
t3 = -g(3) * t16 + t13 * t15;
t35 = pkin(5) * t27 + pkin(7) + pkin(8) + qJ(2);
t12 = g(1) * t28 - g(2) * t30;
t33 = t15 * t17 + t16 * t26;
t25 = cos(pkin(10));
t32 = -t25 * pkin(2) - pkin(1) - t49;
t18 = sin(t23);
t9 = t16 * t37 + t40;
t7 = -t16 * t39 + t38;
t5 = t12 * t15;
t4 = t13 * t16 + t43;
t2 = t3 * t29;
t1 = t3 * t27;
t10 = [0, t12, t13, t12 * t25, -t12 * sin(pkin(10)) -t13, -g(1) * (-t28 * pkin(1) + t30 * qJ(2)) - g(2) * (t30 * pkin(1) + t28 * qJ(2)) 0, 0, 0, 0, 0, t12 * t19, -t12 * t18, 0, 0, 0, 0, 0, t12 * t16, -t5, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t5 (-g(1) * t35 + g(2) * t32) * t30 + (-g(1) * t32 - g(2) * t35) * t28; 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t19 + t13 * t18, g(3) * t18 + t13 * t19, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(3) * t49 + t13 * (pkin(3) * t18 + t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(3) * t34 + t13 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, g(1) * t9 - g(2) * t7 + t29 * t43, 0, t48 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3;];
taug_reg  = t10;
