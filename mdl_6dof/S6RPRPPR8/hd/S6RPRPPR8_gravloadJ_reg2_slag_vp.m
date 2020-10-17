% Calculate inertial parameters regressor of gravitation load for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:23:37
% EndTime: 2019-05-05 17:23:38
% DurationCPUTime: 0.29s
% Computational Cost: add. (133->81), mult. (299->92), div. (0->0), fcn. (288->6), ass. (0->46)
t30 = sin(qJ(1));
t29 = sin(qJ(3));
t33 = cos(qJ(1));
t48 = t29 * t33;
t59 = pkin(4) * t48 + t30 * qJ(5);
t32 = cos(qJ(3));
t22 = t32 * qJ(4);
t58 = t32 * pkin(5) + t22;
t51 = g(3) * t32;
t53 = g(2) * t33;
t54 = g(1) * t30;
t9 = -t53 + t54;
t57 = t9 * t29 + t51;
t56 = -pkin(1) - pkin(7);
t55 = -pkin(3) - pkin(4);
t52 = g(3) * t29;
t49 = t29 * t30;
t47 = t30 * t32;
t28 = sin(qJ(6));
t46 = t33 * t28;
t31 = cos(qJ(6));
t45 = t33 * t31;
t42 = qJ(4) * t29;
t44 = pkin(3) * t47 + t30 * t42;
t43 = t33 * pkin(1) + t30 * qJ(2);
t41 = -pkin(8) + t55;
t40 = pkin(4) * t47 + t44;
t39 = t33 * pkin(7) + t43;
t38 = pkin(3) * t49 + t39;
t23 = t33 * qJ(2);
t37 = t56 * t30 + t23;
t10 = g(1) * t33 + g(2) * t30;
t36 = -t29 * pkin(8) + t58;
t35 = -t30 * t22 + t38;
t17 = pkin(3) * t48;
t34 = -t33 * t22 + t17 + t37;
t12 = pkin(4) * t49;
t8 = t10 * t32;
t7 = t10 * t29;
t6 = t30 * t28 - t32 * t45;
t5 = t30 * t31 + t32 * t46;
t4 = t31 * t47 + t46;
t3 = t28 * t47 - t45;
t2 = g(1) * t47 - t32 * t53 - t52;
t1 = g(1) * t49 - g(2) * t48 + t51;
t11 = [0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10, -g(1) * (-t30 * pkin(1) + t23) - g(2) * t43, 0, 0, 0, 0, 0, 0, -t7, -t8, t9, -g(1) * t37 - g(2) * t39, 0, 0, 0, 0, 0, 0, -t7, t9, t8, -g(1) * t34 - g(2) * t35, 0, 0, 0, 0, 0, 0, t8, t7, -t9, -g(1) * (t34 + t59) - g(2) * (-t33 * qJ(5) + t12 + t35) 0, 0, 0, 0, 0, 0, -g(1) * t6 + g(2) * t4, -g(1) * t5 - g(2) * t3, -t7, -g(1) * (t17 + t23 + t59) - g(2) * (t12 + t38) + (g(1) * t36 + g(2) * qJ(5)) * t33 + (-g(1) * t56 + g(2) * t36) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, -t1, -g(1) * t44 - g(3) * (-t29 * pkin(3) + t22) - (-pkin(3) * t32 - t42) * t53, 0, 0, 0, 0, 0, 0, -t1, t2, 0, -g(1) * t40 - g(3) * (t55 * t29 + t22) - (t55 * t32 - t42) * t53, 0, 0, 0, 0, 0, 0, -t57 * t31, t57 * t28, -t2, -g(1) * (pkin(8) * t47 + t40) - g(3) * t58 + (-pkin(5) * t54 - g(3) * t41) * t29 - ((-pkin(5) - qJ(4)) * t29 + t41 * t32) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t5 + t28 * t52, -g(1) * t4 - g(2) * t6 + t31 * t52, 0, 0;];
taug_reg  = t11;
