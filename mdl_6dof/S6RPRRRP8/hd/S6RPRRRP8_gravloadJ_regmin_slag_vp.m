% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP8
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
% taug_reg [6x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:48:04
% EndTime: 2019-05-06 01:48:05
% DurationCPUTime: 0.26s
% Computational Cost: add. (242->70), mult. (346->91), div. (0->0), fcn. (359->8), ass. (0->46)
t29 = qJ(3) + qJ(4);
t24 = cos(t29);
t21 = t24 * pkin(9);
t31 = sin(qJ(3));
t57 = -pkin(3) * t31 + t21;
t34 = cos(qJ(3));
t55 = pkin(3) * t34;
t35 = cos(qJ(1));
t54 = g(2) * t35;
t23 = sin(t29);
t53 = g(3) * t23;
t52 = g(3) * t24;
t30 = sin(qJ(5));
t51 = g(3) * t30;
t32 = sin(qJ(1));
t50 = t23 * t32;
t49 = t30 * t32;
t48 = t30 * t35;
t33 = cos(qJ(5));
t47 = t32 * t33;
t46 = t35 * t33;
t45 = t35 * pkin(1) + t32 * qJ(2);
t44 = t24 * t49;
t43 = pkin(9) * t50 + qJ(6) * t44 + (pkin(4) * t32 + pkin(5) * t47) * t24;
t10 = t23 * t48 + t47;
t8 = t23 * t49 - t46;
t42 = g(1) * t10 + g(2) * t8;
t16 = g(1) * t35 + g(2) * t32;
t15 = g(1) * t32 - t54;
t41 = -pkin(5) * t33 - qJ(6) * t30 - pkin(4);
t40 = pkin(4) * t23 - t57;
t1 = g(1) * t8 - g(2) * t10 + t24 * t51;
t11 = t23 * t46 - t49;
t9 = t23 * t47 + t48;
t39 = g(1) * t9 - g(2) * t11 + t33 * t52;
t38 = t41 * t53;
t6 = -t15 * t24 + t53;
t37 = -pkin(9) * t23 + t41 * t24;
t36 = -pkin(8) - pkin(7);
t26 = t35 * qJ(2);
t7 = t16 * t24;
t5 = g(1) * t50 - t23 * t54 + t52;
t4 = t6 * t33;
t3 = -g(2) * t24 * t48 + g(1) * t44 - t23 * t51;
t2 = -g(1) * t11 - g(2) * t9;
t12 = [0, t15, t16, -t15, -t16, -g(1) * (-t32 * pkin(1) + t26) - g(2) * t45, 0, 0, 0, 0, 0, -t16 * t31, -t16 * t34, 0, 0, 0, 0, 0, -t16 * t23, -t7, 0, 0, 0, 0, 0, t2, t42, t2, t7, -t42, -g(1) * (t11 * pkin(5) + t10 * qJ(6) + t26) - g(2) * (t9 * pkin(5) + t8 * qJ(6) + t45) + (-g(1) * t40 + g(2) * t36) * t35 + (-g(1) * (-pkin(1) + t36) - g(2) * t40) * t32; 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t31 - t15 * t34, g(3) * t34 + t15 * t31, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t4, t3, t4, -t5, -t3, -g(1) * (t32 * t55 + t43) - g(3) * t57 - t38 - (t37 - t55) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, t4, t3, t4, -t5, -t3, -g(1) * t43 - g(3) * t21 - t37 * t54 - t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t39, t1, 0, -t39, -g(1) * (-pkin(5) * t8 + qJ(6) * t9) - g(2) * (pkin(5) * t10 - qJ(6) * t11) - (-pkin(5) * t30 + qJ(6) * t33) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t12;
