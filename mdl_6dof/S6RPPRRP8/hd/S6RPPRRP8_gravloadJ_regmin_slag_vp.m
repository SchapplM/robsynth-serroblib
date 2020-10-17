% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:07:24
% EndTime: 2019-05-05 15:07:25
% DurationCPUTime: 0.23s
% Computational Cost: add. (186->62), mult. (271->74), div. (0->0), fcn. (285->8), ass. (0->34)
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t43 = -g(1) * t25 + g(2) * t27;
t20 = pkin(9) + qJ(4);
t14 = sin(t20);
t15 = cos(t20);
t42 = -g(3) * t14 - t43 * t15;
t38 = g(3) * t15;
t24 = sin(qJ(5));
t37 = t25 * t24;
t26 = cos(qJ(5));
t36 = t25 * t26;
t35 = t27 * t24;
t34 = t27 * t26;
t33 = t27 * pkin(1) + t25 * qJ(2);
t7 = t14 * t37 - t34;
t9 = t14 * t35 + t36;
t32 = g(1) * t9 + g(2) * t7;
t12 = g(1) * t27 + g(2) * t25;
t31 = pkin(5) * t26 + qJ(6) * t24 + pkin(4);
t21 = sin(pkin(9));
t30 = pkin(3) * t21 + t14 * pkin(4) - t15 * pkin(8);
t1 = g(1) * t7 - g(2) * t9 + t24 * t38;
t10 = t14 * t34 - t37;
t8 = t14 * t36 + t35;
t29 = g(1) * t8 - g(2) * t10 + t26 * t38;
t23 = -pkin(7) - qJ(3);
t17 = t27 * qJ(2);
t6 = t12 * t15;
t5 = -t14 * t43 + t38;
t4 = t42 * t26;
t3 = t42 * t24;
t2 = -g(1) * t10 - g(2) * t8;
t11 = [0, -t43, t12, t43, -t12, -g(1) * (-t25 * pkin(1) + t17) - g(2) * t33, -t12 * t21, -t12 * cos(pkin(9)) -t43, -g(1) * (t17 + (-pkin(1) - qJ(3)) * t25) - g(2) * (t27 * qJ(3) + t33) 0, 0, 0, 0, 0, -t12 * t14, -t6, 0, 0, 0, 0, 0, t2, t32, t2, t6, -t32, -g(1) * (t10 * pkin(5) + t9 * qJ(6) + t17) - g(2) * (t8 * pkin(5) + t7 * qJ(6) + t33) + (-g(1) * t30 + g(2) * t23) * t27 + (-g(1) * (-pkin(1) + t23) - g(2) * t30) * t25; 0, 0, 0, 0, 0, t43, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t5, 0, 0, 0, 0, 0, -t4, t3, -t4, -t5, -t3 (pkin(8) * t43 + g(3) * t31) * t14 + (-g(3) * pkin(8) + t43 * t31) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t29, t1, 0, -t29, -g(1) * (-t7 * pkin(5) + t8 * qJ(6)) - g(2) * (t9 * pkin(5) - t10 * qJ(6)) - (-pkin(5) * t24 + qJ(6) * t26) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t11;
