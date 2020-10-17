% Calculate minimal parameter regressor of gravitation load for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:59:06
% EndTime: 2019-05-04 22:59:06
% DurationCPUTime: 0.24s
% Computational Cost: add. (281->69), mult. (498->108), div. (0->0), fcn. (604->12), ass. (0->40)
t26 = sin(pkin(10));
t31 = sin(qJ(2));
t33 = cos(qJ(2));
t42 = cos(pkin(10));
t43 = cos(pkin(6));
t38 = t43 * t42;
t12 = t26 * t31 - t33 * t38;
t41 = t26 * t43;
t14 = t42 * t31 + t33 * t41;
t27 = sin(pkin(6));
t50 = g(3) * t27;
t4 = -g(1) * t14 - g(2) * t12 + t33 * t50;
t24 = pkin(11) + qJ(4);
t22 = sin(t24);
t30 = sin(qJ(6));
t49 = t22 * t30;
t32 = cos(qJ(6));
t48 = t22 * t32;
t47 = t26 * t27;
t46 = t27 * t31;
t45 = t30 * t33;
t44 = t32 * t33;
t40 = t27 * t42;
t13 = t26 * t33 + t31 * t38;
t15 = -t31 * t41 + t42 * t33;
t39 = g(1) * t15 + g(2) * t13;
t23 = cos(t24);
t10 = t22 * t46 - t43 * t23;
t6 = t13 * t22 + t23 * t40;
t8 = t15 * t22 - t23 * t47;
t36 = g(1) * t8 + g(2) * t6 + g(3) * t10;
t11 = t43 * t22 + t23 * t46;
t7 = t13 * t23 - t22 * t40;
t9 = t15 * t23 + t22 * t47;
t35 = g(1) * t9 + g(2) * t7 + g(3) * t11;
t34 = g(3) * t46 + t39;
t28 = cos(pkin(11));
t3 = t4 * t23;
t2 = t4 * t22;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t4, t34, -t4 * t28, t4 * sin(pkin(11)) -t34, -g(1) * (-t14 * pkin(2) + t15 * qJ(3)) - g(2) * (-t12 * pkin(2) + t13 * qJ(3)) - (pkin(2) * t33 + qJ(3) * t31) * t50, 0, 0, 0, 0, 0, -t3, t2, -t34, t3, -t2 (t31 * t50 + t39) * (-pkin(8) - qJ(3)) - t4 * (t28 * pkin(3) + pkin(4) * t23 + qJ(5) * t22 + pkin(2)) 0, 0, 0, 0, 0, -g(1) * (-t14 * t49 + t15 * t32) - g(2) * (-t12 * t49 + t13 * t32) - (t22 * t45 + t31 * t32) * t50, -g(1) * (-t14 * t48 - t15 * t30) - g(2) * (-t12 * t48 - t13 * t30) - (t22 * t44 - t30 * t31) * t50; 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t35, 0, -t36, -t35, -g(1) * (-t8 * pkin(4) + t9 * qJ(5)) - g(2) * (-t6 * pkin(4) + t7 * qJ(5)) - g(3) * (-t10 * pkin(4) + t11 * qJ(5)) 0, 0, 0, 0, 0, -t35 * t30, -t35 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t30 + t8 * t32) - g(2) * (-t12 * t30 + t6 * t32) - g(3) * (t10 * t32 + t27 * t45) -g(1) * (-t14 * t32 - t8 * t30) - g(2) * (-t12 * t32 - t6 * t30) - g(3) * (-t10 * t30 + t27 * t44);];
taug_reg  = t1;
