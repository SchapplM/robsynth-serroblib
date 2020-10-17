% Calculate inertial parameters regressor of gravitation load for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:21:40
% EndTime: 2019-05-04 23:21:41
% DurationCPUTime: 0.40s
% Computational Cost: add. (290->97), mult. (740->132), div. (0->0), fcn. (890->10), ass. (0->55)
t33 = sin(pkin(10));
t35 = cos(pkin(10));
t41 = cos(qJ(2));
t38 = sin(qJ(2));
t58 = cos(pkin(6));
t54 = t38 * t58;
t21 = t33 * t41 + t35 * t54;
t23 = -t33 * t54 + t35 * t41;
t73 = -g(1) * t23 - g(2) * t21;
t37 = sin(qJ(4));
t72 = pkin(4) * t37;
t34 = sin(pkin(6));
t69 = g(3) * t34;
t53 = t41 * t58;
t20 = t33 * t38 - t35 * t53;
t68 = t20 * pkin(8);
t22 = t33 * t53 + t35 * t38;
t67 = t22 * pkin(8);
t66 = t34 * t37;
t65 = t34 * t38;
t40 = cos(qJ(4));
t64 = t34 * t40;
t63 = t34 * t41;
t36 = sin(qJ(6));
t62 = t36 * t40;
t39 = cos(qJ(6));
t61 = t39 * t40;
t60 = pkin(2) * t63 + qJ(3) * t65;
t59 = qJ(5) * t40;
t57 = pkin(8) * t63 + t60;
t10 = t22 * t37 + t33 * t64;
t9 = -t22 * t40 + t33 * t66;
t56 = -t9 * pkin(4) + qJ(5) * t10;
t11 = t20 * t40 + t35 * t66;
t12 = -t20 * t37 + t35 * t64;
t55 = t11 * pkin(4) - qJ(5) * t12;
t17 = t20 * pkin(2);
t52 = t21 * qJ(3) - t17;
t18 = t22 * pkin(2);
t51 = t23 * qJ(3) - t18;
t24 = t58 * t37 + t40 * t63;
t25 = -t37 * t63 + t58 * t40;
t50 = -t24 * pkin(4) + qJ(5) * t25;
t49 = t65 * t72 + t57;
t48 = qJ(3) - t59;
t47 = t21 * t72 - t17 - t68;
t46 = t23 * t72 - t18 - t67;
t45 = pkin(9) * t37 - t59;
t2 = g(1) * t9 - g(2) * t11 + g(3) * t24;
t43 = g(1) * t10 - g(2) * t12 + g(3) * t25;
t5 = -g(1) * t22 - g(2) * t20 + g(3) * t63;
t42 = g(3) * t65 - t73;
t4 = t42 * t40;
t3 = t42 * t37;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t42, -g(1) * t51 - g(2) * t52 - g(3) * t60, 0, 0, 0, 0, 0, 0, -t3, -t4, -t5, -g(1) * (t51 - t67) - g(2) * (t52 - t68) - g(3) * t57, 0, 0, 0, 0, 0, 0, -t5, t3, t4, -g(1) * (t48 * t23 + t46) - g(2) * (t48 * t21 + t47) - g(3) * (-t59 * t65 + t49) 0, 0, 0, 0, 0, 0, -g(1) * (-t22 * t39 - t23 * t62) - g(2) * (-t20 * t39 - t21 * t62) - (-t38 * t62 + t39 * t41) * t69, -g(1) * (t22 * t36 - t23 * t61) - g(2) * (t20 * t36 - t21 * t61) - (-t36 * t41 - t38 * t61) * t69, -t3, -g(1) * (-t22 * pkin(5) + t46) - g(2) * (-t20 * pkin(5) + t47) - g(3) * t49 - (pkin(5) * t41 + t45 * t38) * t69 + t73 * (qJ(3) + t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t43, -g(1) * t56 - g(2) * t55 - g(3) * t50, 0, 0, 0, 0, 0, 0, -t43 * t36, -t43 * t39, t2, -g(1) * (-pkin(9) * t9 + t56) - g(2) * (pkin(9) * t11 + t55) - g(3) * (-pkin(9) * t24 + t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t23 * t36 + t39 * t9) - g(2) * (-t11 * t39 - t21 * t36) - g(3) * (t24 * t39 - t36 * t65) -g(1) * (-t23 * t39 - t36 * t9) - g(2) * (t11 * t36 - t21 * t39) - g(3) * (-t24 * t36 - t39 * t65) 0, 0;];
taug_reg  = t1;
