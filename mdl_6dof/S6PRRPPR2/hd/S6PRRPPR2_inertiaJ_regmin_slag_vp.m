% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:46:33
% EndTime: 2019-05-05 02:46:34
% DurationCPUTime: 0.42s
% Computational Cost: add. (365->68), mult. (781->134), div. (0->0), fcn. (956->10), ass. (0->55)
t34 = sin(pkin(11));
t36 = cos(pkin(11));
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t21 = t34 * t38 - t36 * t41;
t22 = t34 * t41 + t36 * t38;
t31 = -t41 * pkin(3) - pkin(2);
t46 = -t22 * qJ(5) + t31;
t12 = t21 * pkin(4) + t46;
t65 = -0.2e1 * t12;
t26 = t34 * pkin(3) + qJ(5);
t64 = 0.2e1 * t26;
t63 = 0.2e1 * t41;
t62 = t26 * t21;
t35 = sin(pkin(6));
t61 = t35 * sin(qJ(2));
t42 = cos(qJ(2));
t60 = t35 * t42;
t37 = sin(qJ(6));
t59 = t37 * t21;
t58 = t37 * t22;
t40 = cos(qJ(6));
t57 = t40 * t21;
t18 = t40 * t22;
t56 = t40 * t37;
t55 = -qJ(4) - pkin(8);
t54 = cos(pkin(6));
t53 = 0.2e1 * t21 * t22;
t24 = t55 * t41;
t51 = t55 * t38;
t13 = -t34 * t24 - t36 * t51;
t15 = -t36 * t24 + t34 * t51;
t52 = t13 ^ 2 + t15 ^ 2;
t30 = -t36 * pkin(3) - pkin(4);
t19 = t54 * t38 + t41 * t61;
t44 = -t38 * t61 + t54 * t41;
t10 = t36 * t19 + t34 * t44;
t8 = t34 * t19 - t36 * t44;
t50 = t35 ^ 2 * t42 ^ 2 + t10 ^ 2 + t8 ^ 2;
t49 = t10 * t15 + t8 * t13;
t48 = -t10 * t21 + t8 * t22;
t25 = -pkin(9) + t30;
t47 = -t22 * t25 + t62;
t45 = 0.2e1 * t13 * t22 - 0.2e1 * t15 * t21;
t33 = t40 ^ 2;
t32 = t37 ^ 2;
t20 = t21 ^ 2;
t7 = -t21 * pkin(5) + t15;
t6 = t22 * pkin(5) + t13;
t5 = (pkin(4) + pkin(9)) * t21 + t46;
t4 = -t37 * t8 + t40 * t60;
t3 = t37 * t60 + t40 * t8;
t2 = t37 * t6 + t40 * t5;
t1 = -t37 * t5 + t40 * t6;
t9 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0; 0, 0, t60, -t61, 0, 0, 0, 0, 0, t41 * t60, -t38 * t60, t48, -t31 * t60 + t49, t48, t21 * t60, t22 * t60, -t12 * t60 + t49, 0, 0, 0, 0, 0, -t10 * t57 + t3 * t22, t10 * t59 + t4 * t22; 0, 1, 0, 0, t38 ^ 2, t38 * t63, 0, 0, 0, pkin(2) * t63, -0.2e1 * pkin(2) * t38, t45, t31 ^ 2 + t52, t45, t21 * t65, t22 * t65, t12 ^ 2 + t52, t32 * t20, 0.2e1 * t20 * t56, t37 * t53, t40 * t53, t22 ^ 2, 0.2e1 * t1 * t22 - 0.2e1 * t7 * t57, -0.2e1 * t2 * t22 + 0.2e1 * t7 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t19, 0 (t10 * t34 - t36 * t8) * pkin(3), 0, t8, t10, t10 * t26 + t8 * t30, 0, 0, 0, 0, 0, t10 * t37, t10 * t40; 0, 0, 0, 0, 0, 0, t38, t41, 0, -t38 * pkin(8), -t41 * pkin(8) (-t21 * t34 - t22 * t36) * pkin(3) (-t13 * t36 + t15 * t34) * pkin(3), t30 * t22 - t62, t13, t15, t13 * t30 + t15 * t26, t21 * t56 (-t32 + t33) * t21, t18, -t58, 0, t7 * t37 - t47 * t40, t47 * t37 + t7 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t34 ^ 2 + t36 ^ 2) * pkin(3) ^ 2, 0, 0.2e1 * t30, t64, t26 ^ 2 + t30 ^ 2, t33, -0.2e1 * t56, 0, 0, 0, t37 * t64, t40 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, 0, 0, 0, -t60, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t21, -t22, t12, 0, 0, 0, 0, 0, -t58, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, t13, 0, 0, 0, 0, 0, t18, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t57, t22, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t37, 0, t40 * t25, -t37 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
