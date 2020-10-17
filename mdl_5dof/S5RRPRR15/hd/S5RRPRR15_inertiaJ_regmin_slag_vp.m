% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR15_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:08
% EndTime: 2019-12-31 20:43:10
% DurationCPUTime: 0.36s
% Computational Cost: add. (233->69), mult. (472->126), div. (0->0), fcn. (511->6), ass. (0->61)
t35 = sin(qJ(4));
t24 = t35 * pkin(4) + qJ(3);
t64 = 0.2e1 * t24;
t36 = sin(qJ(2));
t63 = -0.2e1 * t36;
t62 = 0.2e1 * t36;
t39 = cos(qJ(2));
t61 = 0.2e1 * t39;
t60 = 0.2e1 * qJ(3);
t40 = -pkin(2) - pkin(7);
t34 = sin(qJ(5));
t59 = t34 * pkin(4);
t58 = t36 * pkin(4);
t37 = cos(qJ(5));
t57 = t37 * pkin(4);
t38 = cos(qJ(4));
t44 = -t36 * qJ(3) - pkin(1);
t14 = t40 * t39 + t44;
t45 = pkin(8) * t39 - t14;
t27 = t36 * pkin(6);
t21 = t36 * pkin(3) + t27;
t54 = t35 * t21;
t5 = -t45 * t38 + t54;
t56 = t37 * t5;
t15 = t34 * t38 + t37 * t35;
t55 = t15 * t36;
t53 = t35 * t36;
t52 = t35 * t39;
t51 = t36 * t39;
t50 = t38 * t35;
t49 = t38 * t39;
t28 = t39 * pkin(6);
t22 = t39 * pkin(3) + t28;
t31 = t36 ^ 2;
t33 = t39 ^ 2;
t48 = t31 + t33;
t47 = t39 * qJ(3);
t46 = -0.2e1 * t51;
t17 = t38 * t21;
t4 = t45 * t35 + t17 + t58;
t1 = -t34 * t5 + t37 * t4;
t43 = -t36 * pkin(2) + t47;
t42 = t36 * t40 + t47;
t32 = t38 ^ 2;
t30 = t35 ^ 2;
t26 = t38 * t40;
t25 = t38 * t36;
t20 = -t39 * pkin(2) + t44;
t19 = -t38 * pkin(8) + t26;
t18 = (-pkin(8) + t40) * t35;
t16 = -t34 * t35 + t37 * t38;
t13 = t16 * t36;
t12 = pkin(4) * t49 + t22;
t11 = t15 * t39;
t10 = t34 * t52 - t37 * t49;
t9 = t37 * t18 + t34 * t19;
t8 = -t34 * t18 + t37 * t19;
t7 = t38 * t14 + t54;
t6 = -t35 * t14 + t17;
t2 = t34 * t4 + t56;
t3 = [1, 0, 0, t31, 0.2e1 * t51, 0, 0, 0, pkin(1) * t61, pkin(1) * t63, 0.2e1 * t48 * pkin(6), t20 * t61, t20 * t63, t48 * pkin(6) ^ 2 + t20 ^ 2, t30 * t33, 0.2e1 * t33 * t50, t35 * t46, t38 * t46, t31, 0.2e1 * t22 * t49 + 0.2e1 * t6 * t36, -0.2e1 * t22 * t52 - 0.2e1 * t7 * t36, t11 ^ 2, -0.2e1 * t11 * t10, -t11 * t62, t10 * t62, t31, 0.2e1 * t1 * t36 - 0.2e1 * t12 * t10, -0.2e1 * t12 * t11 - 0.2e1 * t2 * t36; 0, 0, 0, 0, 0, t36, t39, 0, -t27, -t28, t43, t27, t28, t43 * pkin(6), -t35 * t49, (t30 - t32) * t39, t25, -t53, 0, t22 * t35 + t42 * t38, t22 * t38 - t42 * t35, -t11 * t16, t16 * t10 + t11 * t15, t13, -t55, 0, -t24 * t10 + t12 * t15 + t8 * t36, -t24 * t11 + t12 * t16 - t9 * t36; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t60, pkin(2) ^ 2 + qJ(3) ^ 2, t32, -0.2e1 * t50, 0, 0, 0, t35 * t60, t38 * t60, t16 ^ 2, -0.2e1 * t16 * t15, 0, 0, 0, t15 * t64, t16 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, t27, 0, 0, 0, 0, 0, t25, -t53, 0, 0, 0, 0, 0, t13, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t49, t36, t6, -t7, 0, 0, -t11, t10, t36, t36 * t57 + t1, -t56 + (-t4 - t58) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, t26, -t35 * t40, 0, 0, t16, -t15, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t35, 0, 0, 0, 0, 0, t16, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, t36, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t57, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
