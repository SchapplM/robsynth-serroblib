% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:06:01
% EndTime: 2019-05-05 16:06:02
% DurationCPUTime: 0.29s
% Computational Cost: add. (267->76), mult. (269->83), div. (0->0), fcn. (258->10), ass. (0->46)
t30 = pkin(10) + qJ(4);
t24 = qJ(5) + t30;
t20 = sin(t24);
t21 = cos(t24);
t43 = -t20 * pkin(5) + t21 * pkin(9);
t35 = sin(qJ(1));
t37 = cos(qJ(1));
t56 = g(2) * t37;
t16 = g(1) * t35 - t56;
t60 = -g(3) * t20 + t16 * t21;
t22 = sin(t30);
t59 = pkin(4) * t22;
t23 = cos(t30);
t58 = pkin(4) * t23;
t57 = pkin(5) * t21;
t54 = g(3) * t21;
t31 = sin(pkin(10));
t52 = t31 * pkin(3);
t51 = t20 * t35;
t34 = sin(qJ(6));
t50 = t35 * t34;
t36 = cos(qJ(6));
t49 = t35 * t36;
t48 = t37 * t34;
t47 = t37 * t36;
t33 = -pkin(7) - qJ(3);
t46 = pkin(9) * t51 + t35 * t57;
t45 = t37 * pkin(1) + t35 * qJ(2);
t26 = t37 * qJ(2);
t44 = -t35 * pkin(1) + t26;
t42 = -pkin(9) * t20 - t57;
t17 = g(1) * t37 + g(2) * t35;
t12 = t52 + t59;
t29 = -pkin(8) + t33;
t40 = t37 * t12 + t35 * t29 + t44;
t39 = t35 * t12 - t37 * t29 + t45;
t38 = g(3) * t22 - t16 * t23;
t9 = t20 * t47 - t50;
t8 = t20 * t48 + t49;
t7 = t20 * t49 + t48;
t6 = -t20 * t50 + t47;
t5 = t17 * t21;
t3 = g(1) * t51 - t20 * t56 + t54;
t2 = t60 * t36;
t1 = t60 * t34;
t4 = [0, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, -g(1) * t44 - g(2) * t45, 0, 0, 0, 0, 0, 0, -t17 * t31, -t17 * cos(pkin(10)) t16, -g(1) * (t26 + (-pkin(1) - qJ(3)) * t35) - g(2) * (t37 * qJ(3) + t45) 0, 0, 0, 0, 0, 0, -t17 * t22, -t17 * t23, t16, -g(1) * (t37 * t52 + t26 + (-pkin(1) + t33) * t35) - g(2) * (-t37 * t33 + t35 * t52 + t45) 0, 0, 0, 0, 0, 0, -t17 * t20, -t5, t16, -g(1) * t40 - g(2) * t39, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7, g(1) * t8 - g(2) * t6, t5, -g(1) * (-t37 * t43 + t40) - g(2) * (-t35 * t43 + t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, g(3) * t23 + t16 * t22, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t3, 0, t38 * pkin(4), 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * (t35 * t58 + t46) - g(3) * (t43 - t59) - (t42 - t58) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * t46 - g(3) * t43 - t42 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8 + t34 * t54, g(1) * t7 - g(2) * t9 + t36 * t54, 0, 0;];
taug_reg  = t4;
