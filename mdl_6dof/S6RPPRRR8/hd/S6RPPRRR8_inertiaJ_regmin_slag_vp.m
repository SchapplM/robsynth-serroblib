% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:15:02
% EndTime: 2019-05-05 16:15:04
% DurationCPUTime: 0.58s
% Computational Cost: add. (506->83), mult. (967->137), div. (0->0), fcn. (1208->8), ass. (0->64)
t45 = sin(pkin(10));
t46 = cos(pkin(10));
t50 = sin(qJ(4));
t63 = cos(qJ(4));
t26 = t63 * t45 + t50 * t46;
t22 = t26 ^ 2;
t25 = t50 * t45 - t63 * t46;
t23 = t25 ^ 2;
t75 = -t22 - t23;
t48 = sin(qJ(6));
t49 = sin(qJ(5));
t51 = cos(qJ(6));
t52 = cos(qJ(5));
t29 = t48 * t52 + t51 * t49;
t62 = t25 * t29;
t74 = 0.2e1 * t62;
t73 = -0.2e1 * t25;
t72 = 0.2e1 * t26;
t38 = -t52 * pkin(5) - pkin(4);
t71 = 0.2e1 * t38;
t70 = 2 * qJ(2);
t69 = pkin(8) + pkin(9);
t68 = t26 * pkin(5);
t67 = t48 * pkin(5);
t66 = t51 * pkin(5);
t35 = t45 * pkin(3) + qJ(2);
t14 = t26 * pkin(4) + t25 * pkin(8) + t35;
t47 = -pkin(1) - qJ(3);
t64 = -pkin(7) + t47;
t30 = t64 * t45;
t31 = t64 * t46;
t16 = t63 * t30 + t50 * t31;
t56 = t52 * t16;
t5 = t56 + (pkin(9) * t25 + t14) * t49;
t65 = t51 * t5;
t28 = t48 * t49 - t51 * t52;
t13 = t25 * t28;
t61 = t25 * t49;
t60 = t25 * t52;
t59 = t29 * t26;
t58 = t49 * t26;
t57 = t49 * t52;
t21 = t52 * t26;
t34 = t45 ^ 2 + t46 ^ 2;
t55 = t25 * t72;
t6 = t52 * t14 - t49 * t16;
t4 = pkin(9) * t60 + t6 + t68;
t1 = t51 * t4 - t48 * t5;
t54 = pkin(4) * t25 - pkin(8) * t26;
t15 = t50 * t30 - t63 * t31;
t53 = qJ(2) ^ 2;
t44 = t52 ^ 2;
t43 = t49 ^ 2;
t33 = t69 * t52;
t32 = t69 * t49;
t24 = t34 * t47;
t19 = -t48 * t32 + t51 * t33;
t18 = -t51 * t32 - t48 * t33;
t17 = t28 * t26;
t12 = t51 * t21 - t48 * t58;
t8 = -pkin(5) * t61 + t15;
t7 = t49 * t14 + t56;
t2 = t48 * t4 + t65;
t3 = [1, 0, 0, -2 * pkin(1), t70, pkin(1) ^ 2 + t53, t45 * t70, t46 * t70, -0.2e1 * t24, t34 * t47 ^ 2 + t53, t23, t55, 0, 0, 0, t35 * t72, t35 * t73, t44 * t23, -0.2e1 * t23 * t57, t21 * t73, t49 * t55, t22, -0.2e1 * t15 * t61 + 0.2e1 * t6 * t26, -0.2e1 * t15 * t60 - 0.2e1 * t7 * t26, t13 ^ 2, t13 * t74, t13 * t72, t26 * t74, t22, 0.2e1 * t1 * t26 - 0.2e1 * t62 * t8, 0.2e1 * t8 * t13 - 0.2e1 * t2 * t26; 0, 0, 0, 1, 0, -pkin(1), 0, 0, -t34, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75 * t49, t75 * t52, 0, 0, 0, 0, 0, -t25 * t62 - t26 * t59, -t12 * t26 + t25 * t13; 0, 0, 0, 0, 0, 1, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t45, t46, 0, qJ(2), 0, 0, 0, 0, 0, t26, -t25, 0, 0, 0, 0, 0, t21, -t58, 0, 0, 0, 0, 0, -t17, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, 0, -t15, -t16, -t25 * t57 (t43 - t44) * t25, t58, t21, 0, -t15 * t52 + t54 * t49, t15 * t49 + t54 * t52, t13 * t29, -t13 * t28 + t29 * t62, t59, -t17, 0, t18 * t26 + t8 * t28 - t38 * t62, t38 * t13 - t19 * t26 + t8 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t26, 0, 0, 0, 0, 0, -t60, t61, 0, 0, 0, 0, 0, t13, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t43, 0.2e1 * t57, 0, 0, 0, 0.2e1 * pkin(4) * t52, -0.2e1 * pkin(4) * t49, t29 ^ 2, -0.2e1 * t29 * t28, 0, 0, 0, t28 * t71, t29 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t61, t26, t6, -t7, 0, 0, t13, t62, t26, t26 * t66 + t1, -t65 + (-t4 - t68) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t21, 0, 0, 0, 0, 0, -t59, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t49, 0, 0, 0, 0, 0, -t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t52, 0, -t49 * pkin(8), -t52 * pkin(8), 0, 0, t29, -t28, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t66, -0.2e1 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t62, t26, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t66, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
