% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:15
% EndTime: 2019-12-31 19:39:16
% DurationCPUTime: 0.29s
% Computational Cost: add. (165->73), mult. (313->104), div. (0->0), fcn. (322->8), ass. (0->51)
t36 = sin(qJ(1));
t38 = cos(qJ(1));
t17 = g(1) * t38 + g(2) * t36;
t35 = sin(qJ(2));
t63 = t17 * t35;
t25 = t35 * qJ(3);
t37 = cos(qJ(2));
t51 = t37 * pkin(2) + t25;
t61 = g(1) * t36;
t58 = t37 * pkin(3);
t32 = sin(pkin(8));
t56 = t35 * t32;
t55 = t35 * t38;
t33 = cos(pkin(8));
t22 = t33 * pkin(4) + pkin(3);
t54 = t37 * t22;
t53 = t37 * t32;
t52 = t37 * t38;
t50 = t38 * pkin(1) + t36 * pkin(6);
t49 = qJ(3) * t37;
t48 = pkin(4) * t53;
t47 = pkin(4) * t56;
t46 = pkin(2) * t52 + t38 * t25 + t50;
t16 = -g(2) * t38 + t61;
t31 = pkin(8) + qJ(5);
t23 = sin(t31);
t24 = cos(t31);
t45 = t37 * t23 - t35 * t24;
t44 = t35 * t23 + t37 * t24;
t43 = -t35 * t33 + t53;
t42 = t37 * t33 + t56;
t41 = -pkin(1) - t51;
t1 = t45 * t36;
t3 = t23 * t52 - t24 * t55;
t40 = g(1) * t3 + g(2) * t1 + g(3) * t44;
t2 = t44 * t36;
t4 = t44 * t38;
t39 = g(1) * t4 + g(2) * t2 - g(3) * t45;
t34 = -pkin(7) - qJ(4);
t28 = t38 * pkin(6);
t20 = t38 * t49;
t18 = t36 * t49;
t14 = t16 * t37;
t13 = t16 * t35;
t10 = t42 * t38;
t9 = t43 * t38;
t8 = t42 * t36;
t7 = t43 * t36;
t6 = g(3) * t35 + t17 * t37;
t5 = -g(3) * t37 + t63;
t11 = [0, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, -t17, -g(1) * (-t36 * pkin(1) + t28) - g(2) * t50, 0, 0, 0, 0, 0, 0, t14, -t17, t13, -g(1) * t28 - g(2) * t46 - t41 * t61, 0, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t10, -g(1) * t7 + g(2) * t9, t17, -g(1) * (-t38 * qJ(4) + t28) - g(2) * (pkin(3) * t52 + t46) + (-g(1) * (t41 - t58) + g(2) * qJ(4)) * t36, 0, 0, 0, 0, 0, 0, g(1) * t2 - g(2) * t4, -g(1) * t1 + g(2) * t3, t17, -g(1) * (t38 * t34 + t28) - g(2) * (t22 * t52 + t38 * t47 + t46) + (-g(1) * (t41 - t47 - t54) - g(2) * t34) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, -t6, -g(1) * (-pkin(2) * t55 + t20) - g(2) * (-t36 * t35 * pkin(2) + t18) - g(3) * t51, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t7 - g(3) * t42, -g(1) * t10 - g(2) * t8 + g(3) * t43, 0, -g(1) * t20 - g(2) * t18 - g(3) * (t51 + t58) + (pkin(2) + pkin(3)) * t63, 0, 0, 0, 0, 0, 0, -t40, -t39, 0, -g(1) * (t38 * t48 + t20) - g(2) * (t36 * t48 + t18) - g(3) * (t51 + t54) + (-g(3) * pkin(4) * t32 + t17 * (pkin(2) + t22)) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t39, 0, 0;];
taug_reg = t11;
