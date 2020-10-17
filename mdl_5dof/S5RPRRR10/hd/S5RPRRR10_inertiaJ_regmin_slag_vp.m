% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:44
% EndTime: 2019-12-31 19:10:46
% DurationCPUTime: 0.37s
% Computational Cost: add. (393->71), mult. (831->129), div. (0->0), fcn. (1014->8), ass. (0->58)
t37 = sin(pkin(9));
t38 = cos(pkin(9));
t41 = sin(qJ(3));
t55 = cos(qJ(3));
t21 = t41 * t37 - t55 * t38;
t64 = -0.2e1 * t21;
t63 = 0.2e1 * t21;
t30 = -t38 * pkin(2) - pkin(1);
t62 = 0.2e1 * t30;
t43 = cos(qJ(4));
t32 = -t43 * pkin(4) - pkin(3);
t61 = 0.2e1 * t32;
t60 = pkin(7) + pkin(8);
t59 = t21 * pkin(4);
t39 = sin(qJ(5));
t58 = t39 * pkin(4);
t42 = cos(qJ(5));
t57 = t42 * pkin(4);
t22 = t55 * t37 + t41 * t38;
t12 = t21 * pkin(3) - t22 * pkin(7) + t30;
t40 = sin(qJ(4));
t48 = pkin(6) + qJ(2);
t25 = t48 * t37;
t26 = t48 * t38;
t14 = -t41 * t25 + t55 * t26;
t50 = t43 * t14;
t5 = t50 + (-pkin(8) * t22 + t12) * t40;
t56 = t42 * t5;
t24 = t39 * t43 + t42 * t40;
t54 = t24 * t21;
t53 = t40 * t21;
t52 = t40 * t22;
t51 = t40 * t43;
t49 = t43 * t22;
t47 = t37 ^ 2 + t38 ^ 2;
t46 = t22 * t64;
t6 = t43 * t12 - t40 * t14;
t4 = -pkin(8) * t49 + t59 + t6;
t1 = -t39 * t5 + t42 * t4;
t45 = -pkin(3) * t22 - pkin(7) * t21;
t23 = t39 * t40 - t42 * t43;
t13 = t55 * t25 + t41 * t26;
t36 = t43 ^ 2;
t35 = t40 ^ 2;
t28 = t60 * t43;
t27 = t60 * t40;
t20 = t22 ^ 2;
t19 = t21 ^ 2;
t18 = t43 * t21;
t17 = -t39 * t27 + t42 * t28;
t16 = -t42 * t27 - t39 * t28;
t15 = t23 * t21;
t10 = t23 * t22;
t9 = t24 * t22;
t8 = pkin(4) * t52 + t13;
t7 = t40 * t12 + t50;
t2 = t39 * t4 + t56;
t3 = [1, 0, 0, 0.2e1 * pkin(1) * t38, -0.2e1 * pkin(1) * t37, 0.2e1 * t47 * qJ(2), t47 * qJ(2) ^ 2 + pkin(1) ^ 2, t20, t46, 0, 0, 0, t21 * t62, t22 * t62, t36 * t20, -0.2e1 * t20 * t51, t49 * t63, t40 * t46, t19, 0.2e1 * t13 * t52 + 0.2e1 * t6 * t21, 0.2e1 * t13 * t49 - 0.2e1 * t7 * t21, t10 ^ 2, 0.2e1 * t10 * t9, -t10 * t63, t9 * t64, t19, 0.2e1 * t1 * t21 + 0.2e1 * t8 * t9, -0.2e1 * t8 * t10 - 0.2e1 * t2 * t21; 0, 0, 0, -t38, t37, 0, -pkin(1), 0, 0, 0, 0, 0, t21, t22, 0, 0, 0, 0, 0, t18, -t53, 0, 0, 0, 0, 0, -t15, -t54; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, 0, -t13, -t14, t40 * t49, (-t35 + t36) * t22, t53, t18, 0, -t13 * t43 + t45 * t40, t13 * t40 + t45 * t43, -t10 * t24, t10 * t23 - t24 * t9, t54, -t15, 0, t16 * t21 + t8 * t23 + t32 * t9, -t32 * t10 - t17 * t21 + t8 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t35, 0.2e1 * t51, 0, 0, 0, 0.2e1 * pkin(3) * t43, -0.2e1 * pkin(3) * t40, t24 ^ 2, -0.2e1 * t24 * t23, 0, 0, 0, t23 * t61, t24 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t52, t21, t6, -t7, 0, 0, -t10, -t9, t21, t21 * t57 + t1, -t56 + (-t4 - t59) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t40, 0, 0, 0, 0, 0, -t23, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t40 * pkin(7), -t43 * pkin(7), 0, 0, t24, -t23, 0, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t9, t21, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t23, 0, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t57, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
