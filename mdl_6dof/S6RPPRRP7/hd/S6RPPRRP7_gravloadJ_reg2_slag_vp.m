% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRP7
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:04:26
% EndTime: 2019-05-05 15:04:27
% DurationCPUTime: 0.29s
% Computational Cost: add. (214->64), mult. (299->76), div. (0->0), fcn. (296->8), ass. (0->43)
t31 = -pkin(7) - qJ(3);
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t28 = sin(pkin(9));
t57 = pkin(3) * t28;
t59 = t33 * t31 + t35 * t57;
t58 = -g(1) * t33 + g(2) * t35;
t27 = pkin(9) + qJ(4);
t21 = sin(t27);
t34 = cos(qJ(5));
t46 = t35 * t34;
t32 = sin(qJ(5));
t49 = t33 * t32;
t10 = -t21 * t49 + t46;
t47 = t35 * t32;
t48 = t33 * t34;
t12 = t21 * t47 + t48;
t22 = cos(t27);
t50 = g(3) * t22;
t1 = -g(1) * t10 - g(2) * t12 + t32 * t50;
t8 = -g(3) * t21 - t22 * t58;
t56 = pkin(5) * t32;
t45 = t35 * pkin(1) + t33 * qJ(2);
t43 = t33 * t57 + t45;
t24 = t35 * qJ(2);
t42 = -t33 * pkin(1) + t24;
t40 = t21 * pkin(4) - t22 * pkin(8);
t15 = g(1) * t35 + g(2) * t33;
t20 = t34 * pkin(5) + pkin(4);
t30 = -qJ(6) - pkin(8);
t38 = t21 * t20 + t22 * t30;
t37 = t42 + t59;
t36 = -t35 * t31 + t43;
t13 = t21 * t46 - t49;
t11 = t21 * t48 + t47;
t9 = t15 * t22;
t7 = -t58 * t21 + t50;
t6 = t8 * t34;
t5 = t8 * t32;
t4 = -g(1) * t13 - g(2) * t11;
t3 = g(1) * t12 - g(2) * t10;
t2 = g(1) * t11 - g(2) * t13 + t34 * t50;
t14 = [0, 0, 0, 0, 0, 0, -t58, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t15, -g(1) * t42 - g(2) * t45, 0, 0, 0, 0, 0, 0, -t15 * t28, -t15 * cos(pkin(9)) -t58, -g(1) * (t24 + (-pkin(1) - qJ(3)) * t33) - g(2) * (t35 * qJ(3) + t45) 0, 0, 0, 0, 0, 0, -t15 * t21, -t9, -t58, -g(1) * t37 - g(2) * t36, 0, 0, 0, 0, 0, 0, t4, t3, t9, -g(1) * (t40 * t35 + t37) - g(2) * (t40 * t33 + t36) 0, 0, 0, 0, 0, 0, t4, t3, t9, -g(1) * (t24 + t59) - g(2) * t43 + (-g(1) * t38 - g(2) * (-t31 + t56)) * t35 + (-g(1) * (-pkin(1) - t56) - g(2) * t38) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, -t7, g(3) * t40 + t58 * (pkin(4) * t22 + pkin(8) * t21) 0, 0, 0, 0, 0, 0, -t6, t5, -t7, g(3) * t38 + t58 * (t20 * t22 - t21 * t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg  = t14;
