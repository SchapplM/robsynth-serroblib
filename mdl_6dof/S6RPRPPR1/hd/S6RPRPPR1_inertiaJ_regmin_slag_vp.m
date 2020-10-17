% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:25:22
% EndTime: 2019-05-05 16:25:23
% DurationCPUTime: 0.42s
% Computational Cost: add. (659->80), mult. (1224->153), div. (0->0), fcn. (1431->10), ass. (0->60)
t46 = sin(pkin(9));
t39 = t46 * pkin(1) + pkin(7);
t67 = cos(qJ(3));
t60 = t67 * t39;
t25 = t67 * qJ(4) + t60;
t45 = sin(pkin(10));
t48 = cos(pkin(10));
t51 = sin(qJ(3));
t58 = (-qJ(4) - t39) * t51;
t16 = t45 * t25 - t48 * t58;
t72 = t16 ^ 2;
t28 = t45 * t51 - t48 * t67;
t26 = t28 ^ 2;
t31 = t45 * t67 + t48 * t51;
t44 = sin(pkin(11));
t47 = cos(pkin(11));
t50 = sin(qJ(6));
t52 = cos(qJ(6));
t32 = t52 * t44 + t50 * t47;
t12 = t32 * t31;
t71 = -0.2e1 * t12;
t40 = -t48 * pkin(3) - pkin(4);
t33 = -t47 * pkin(5) + t40;
t70 = 0.2e1 * t33;
t69 = 0.2e1 * t51;
t36 = t45 * pkin(3) + qJ(5);
t68 = pkin(8) + t36;
t49 = cos(pkin(9));
t41 = -t49 * pkin(1) - pkin(2);
t34 = -t67 * pkin(3) + t41;
t11 = t28 * pkin(4) - t31 * qJ(5) + t34;
t18 = t48 * t25 + t45 * t58;
t6 = t44 * t11 + t47 * t18;
t66 = t16 * t28;
t30 = t50 * t44 - t52 * t47;
t20 = t28 * t30;
t65 = t28 * t44;
t64 = t28 * t47;
t21 = t32 * t28;
t63 = t44 * t31;
t62 = t47 * t31;
t61 = t44 ^ 2 + t47 ^ 2;
t5 = t47 * t11 - t44 * t18;
t59 = t61 * t36;
t57 = t6 * t44 + t5 * t47;
t56 = -t5 * t44 + t6 * t47;
t55 = -t28 * t36 + t31 * t40;
t27 = t31 ^ 2;
t24 = t68 * t47;
t23 = t68 * t44;
t19 = t61 * t31;
t15 = -t50 * t23 + t52 * t24;
t14 = -t52 * t23 - t50 * t24;
t13 = t30 * t31;
t7 = pkin(5) * t63 + t16;
t4 = -pkin(8) * t63 + t6;
t3 = t28 * pkin(5) - pkin(8) * t62 + t5;
t2 = t50 * t3 + t52 * t4;
t1 = t52 * t3 - t50 * t4;
t8 = [1, 0, 0 (t46 ^ 2 + t49 ^ 2) * pkin(1) ^ 2, t51 ^ 2, t67 * t69, 0, 0, 0, -0.2e1 * t41 * t67, t41 * t69, 0.2e1 * t16 * t31 - 0.2e1 * t18 * t28, t18 ^ 2 + t34 ^ 2 + t72, 0.2e1 * t16 * t63 + 0.2e1 * t5 * t28, 0.2e1 * t16 * t62 - 0.2e1 * t6 * t28, -0.2e1 * t57 * t31, t5 ^ 2 + t6 ^ 2 + t72, t13 ^ 2, -t13 * t71, -0.2e1 * t13 * t28, t28 * t71, t26, 0.2e1 * t1 * t28 + 0.2e1 * t7 * t12, -0.2e1 * t7 * t13 - 0.2e1 * t2 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t31 + t66, 0, 0, 0, t56 * t31 + t66, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t27 + t26, 0, 0, 0, t61 * t27 + t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t51, t67, 0, -t51 * t39, -t60 (-t28 * t45 - t31 * t48) * pkin(3) (-t16 * t48 + t18 * t45) * pkin(3), -t16 * t47 + t55 * t44, t16 * t44 + t55 * t47, t56, t16 * t40 + t56 * t36, -t13 * t32, -t32 * t12 + t13 * t30, t21, -t20, 0, t33 * t12 + t14 * t28 + t7 * t30, -t33 * t13 - t15 * t28 + t7 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t51, 0 (-t28 * t48 + t31 * t45) * pkin(3), -t64, t65, t19, t28 * t40 + t31 * t59, 0, 0, 0, 0, 0, t20, t21; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t45 ^ 2 + t48 ^ 2) * pkin(3) ^ 2, -0.2e1 * t40 * t47, 0.2e1 * t40 * t44, 0.2e1 * t59, t61 * t36 ^ 2 + t40 ^ 2, t32 ^ 2, -0.2e1 * t32 * t30, 0, 0, 0, t30 * t70, t32 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t64, -t65, -t19, t57, 0, 0, 0, 0, 0, -t20, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t62, 0, t16, 0, 0, 0, 0, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t44, 0, t40, 0, 0, 0, 0, 0, t30, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t12, t28, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t30, 0, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
