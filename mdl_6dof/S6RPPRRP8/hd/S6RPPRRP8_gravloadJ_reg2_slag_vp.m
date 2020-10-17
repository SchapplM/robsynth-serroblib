% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:07:24
% EndTime: 2019-05-05 15:07:26
% DurationCPUTime: 0.32s
% Computational Cost: add. (236->80), mult. (335->91), div. (0->0), fcn. (342->8), ass. (0->48)
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t58 = g(2) * t36;
t11 = g(1) * t34 - t58;
t29 = pkin(9) + qJ(4);
t24 = cos(t29);
t23 = sin(t29);
t57 = g(3) * t23;
t61 = t11 * t24 - t57;
t30 = sin(pkin(9));
t60 = pkin(3) * t30;
t59 = pkin(8) * t23;
t56 = g(3) * t24;
t20 = t24 * pkin(8);
t55 = t23 * t34;
t54 = t23 * t36;
t53 = t24 * t34;
t33 = sin(qJ(5));
t52 = t33 * t36;
t51 = t34 * t33;
t35 = cos(qJ(5));
t50 = t34 * t35;
t49 = t36 * t35;
t48 = pkin(4) * t53 + pkin(8) * t55;
t47 = t36 * pkin(1) + t34 * qJ(2);
t26 = t36 * qJ(2);
t46 = -t34 * pkin(1) + t26;
t7 = t23 * t51 - t49;
t9 = t23 * t52 + t50;
t45 = g(1) * t9 + g(2) * t7;
t12 = g(1) * t36 + g(2) * t34;
t44 = pkin(5) * t35 + qJ(6) * t33;
t32 = -pkin(7) - qJ(3);
t43 = t34 * t32 + t36 * t60 + t46;
t42 = -t32 * t36 + t34 * t60 + t47;
t41 = -pkin(4) - t44;
t1 = g(1) * t7 - g(2) * t9 + t33 * t56;
t10 = t23 * t49 - t51;
t8 = t23 * t50 + t52;
t40 = g(1) * t8 - g(2) * t10 + t35 * t56;
t38 = pkin(4) * t54 - t36 * t20 + t43;
t37 = pkin(4) * t55 - pkin(8) * t53 + t42;
t6 = t12 * t24;
t5 = g(1) * t55 - g(2) * t54 + t56;
t4 = t61 * t35;
t3 = t61 * t33;
t2 = -g(1) * t10 - g(2) * t8;
t13 = [0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -g(1) * t46 - g(2) * t47, 0, 0, 0, 0, 0, 0, -t12 * t30, -t12 * cos(pkin(9)) t11, -g(1) * (t26 + (-pkin(1) - qJ(3)) * t34) - g(2) * (qJ(3) * t36 + t47) 0, 0, 0, 0, 0, 0, -t12 * t23, -t6, t11, -g(1) * t43 - g(2) * t42, 0, 0, 0, 0, 0, 0, t2, t45, t6, -g(1) * t38 - g(2) * t37, 0, 0, 0, 0, 0, 0, t2, t6, -t45, -g(1) * (t10 * pkin(5) + t9 * qJ(6) + t38) - g(2) * (t8 * pkin(5) + t7 * qJ(6) + t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(1) * t48 - g(3) * (-pkin(4) * t23 + t20) - (-pkin(4) * t24 - t59) * t58, 0, 0, 0, 0, 0, 0, -t4, -t5, -t3, -g(1) * (t44 * t53 + t48) - g(3) * t20 - t41 * t57 - (t41 * t24 - t59) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t40, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t40, -g(1) * (-pkin(5) * t7 + qJ(6) * t8) - g(2) * (pkin(5) * t9 - qJ(6) * t10) - (-pkin(5) * t33 + qJ(6) * t35) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;
