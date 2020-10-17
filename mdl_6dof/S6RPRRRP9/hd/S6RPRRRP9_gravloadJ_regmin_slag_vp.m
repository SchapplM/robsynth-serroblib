% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:53:24
% EndTime: 2019-05-06 01:53:24
% DurationCPUTime: 0.22s
% Computational Cost: add. (166->54), mult. (254->82), div. (0->0), fcn. (261->8), ass. (0->43)
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t56 = -g(1) * t31 + g(2) * t34;
t28 = qJ(4) + qJ(5);
t20 = sin(t28);
t30 = sin(qJ(3));
t21 = cos(t28);
t42 = t34 * t21;
t47 = t31 * t20;
t3 = -t30 * t47 + t42;
t43 = t34 * t20;
t46 = t31 * t21;
t5 = t30 * t43 + t46;
t33 = cos(qJ(3));
t50 = g(3) * t33;
t1 = -g(1) * t3 - g(2) * t5 + t20 * t50;
t8 = -g(3) * t30 - t33 * t56;
t29 = sin(qJ(4));
t15 = t29 * pkin(4) + pkin(5) * t20;
t49 = pkin(7) + t15;
t48 = t15 * t30;
t45 = t31 * t29;
t32 = cos(qJ(4));
t44 = t31 * t32;
t41 = t34 * t29;
t40 = t34 * t32;
t16 = t32 * pkin(4) + pkin(5) * t21;
t38 = g(2) * (t34 * pkin(1) + t31 * qJ(2));
t18 = g(1) * t34 + g(2) * t31;
t14 = pkin(3) + t16;
t27 = -qJ(6) - pkin(9) - pkin(8);
t36 = t30 * t14 + t33 * t27;
t23 = t34 * qJ(2);
t13 = t18 * t33;
t12 = t30 * t40 - t45;
t11 = t30 * t41 + t44;
t10 = t30 * t44 + t41;
t9 = -t30 * t45 + t40;
t7 = -t30 * t56 + t50;
t6 = t30 * t42 - t47;
t4 = t30 * t46 + t43;
t2 = g(1) * t4 - g(2) * t6 + t21 * t50;
t17 = [0, -t56, t18, t56, -t18, -g(1) * (-t31 * pkin(1) + t23) - t38, 0, 0, 0, 0, 0, -t18 * t30, -t13, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t13, -g(1) * t23 - t38 + (-g(1) * t36 - g(2) * t49) * t34 + (-g(1) * (-pkin(1) - t49) - g(2) * t36) * t31; 0, 0, 0, 0, 0, t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, -t8 * t32, t8 * t29, 0, 0, 0, 0, 0, -t8 * t21, t8 * t20, -t7, g(3) * t36 + t56 * (t14 * t33 - t27 * t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11 + t29 * t50, g(1) * t10 - g(2) * t12 + t32 * t50, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t34 * t16 - t31 * t48) - g(2) * (t31 * t16 + t34 * t48) + t15 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t17;
