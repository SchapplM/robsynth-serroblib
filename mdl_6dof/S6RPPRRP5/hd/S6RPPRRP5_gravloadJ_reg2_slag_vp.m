% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:58:26
% EndTime: 2019-05-05 14:58:27
% DurationCPUTime: 0.27s
% Computational Cost: add. (129->66), mult. (287->75), div. (0->0), fcn. (284->6), ass. (0->42)
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t15 = g(1) * t29 + g(2) * t26;
t25 = sin(qJ(4));
t24 = sin(qJ(5));
t42 = t29 * t24;
t27 = cos(qJ(5));
t44 = t26 * t27;
t11 = -t25 * t42 - t44;
t28 = cos(qJ(4));
t48 = g(3) * t28;
t41 = t29 * t27;
t45 = t26 * t24;
t9 = t25 * t45 - t41;
t1 = -g(1) * t11 + g(2) * t9 + t24 * t48;
t8 = -g(3) * t25 + t15 * t28;
t52 = g(1) * t26;
t47 = t25 * pkin(4);
t46 = t29 * pkin(7);
t43 = t28 * t29;
t40 = -pkin(1) - qJ(3);
t39 = t29 * pkin(1) + t26 * qJ(2);
t37 = t29 * qJ(3) + t39;
t36 = pkin(5) * t24 + pkin(7);
t35 = g(2) * t37;
t33 = t28 * pkin(8) - t47;
t14 = -g(2) * t29 + t52;
t20 = t29 * qJ(2);
t32 = t40 * t26 + t20;
t17 = t27 * pkin(5) + pkin(4);
t23 = -qJ(6) - pkin(8);
t30 = t25 * t17 + t28 * t23;
t13 = -g(2) * t43 + t28 * t52;
t12 = t25 * t41 - t45;
t10 = -t25 * t44 - t42;
t7 = t15 * t25 + t48;
t6 = t8 * t27;
t5 = t8 * t24;
t4 = -g(1) * t10 - g(2) * t12;
t3 = -g(1) * t9 - g(2) * t11;
t2 = g(1) * t12 - g(2) * t10 + t27 * t48;
t16 = [0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, -g(1) * (-t26 * pkin(1) + t20) - g(2) * t39, 0, 0, 0, 0, 0, 0, 0, -t15, t14, -g(1) * t32 - t35, 0, 0, 0, 0, 0, 0, t14 * t25, t13, t15, -g(1) * (t32 - t46) - g(2) * (-t26 * pkin(7) + t37) 0, 0, 0, 0, 0, 0, t4, t3, -t13, -g(1) * (t20 - t46) - g(2) * (-pkin(8) * t43 + t29 * t47 + t37) + (-g(1) * (t33 + t40) + g(2) * pkin(7)) * t26, 0, 0, 0, 0, 0, 0, t4, t3, -t13, -g(1) * t20 - t35 + (g(1) * t36 - g(2) * t30) * t29 + (-g(1) * (-t30 + t40) + g(2) * t36) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(3) * t33 - t15 * (pkin(4) * t28 + pkin(8) * t25) 0, 0, 0, 0, 0, 0, -t6, t5, -t7, g(3) * t30 - t15 * (t17 * t28 - t23 * t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t16;
