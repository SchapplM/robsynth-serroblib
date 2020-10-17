% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:24:30
% EndTime: 2019-05-06 08:24:32
% DurationCPUTime: 0.58s
% Computational Cost: add. (678->82), mult. (1247->151), div. (0->0), fcn. (1463->8), ass. (0->63)
t49 = sin(pkin(9));
t51 = cos(pkin(9));
t53 = sin(qJ(2));
t55 = cos(qJ(2));
t30 = t49 * t53 - t51 * t55;
t48 = sin(pkin(10));
t50 = cos(pkin(10));
t52 = sin(qJ(6));
t54 = cos(qJ(6));
t31 = t54 * t48 + t52 * t50;
t16 = t31 * t30;
t77 = 0.2e1 * t16;
t33 = -t52 * t48 + t54 * t50;
t41 = t49 * pkin(2) + qJ(4);
t76 = t41 ^ 2;
t32 = t49 * t55 + t51 * t53;
t75 = -0.2e1 * t32;
t35 = t48 * pkin(5) + t41;
t74 = 0.2e1 * t35;
t73 = 0.2e1 * t41;
t72 = 0.2e1 * t55;
t44 = -t51 * pkin(2) - pkin(3);
t40 = -qJ(5) + t44;
t71 = -pkin(8) + t40;
t45 = -t55 * pkin(2) - pkin(1);
t58 = -t32 * qJ(4) + t45;
t10 = (pkin(3) + qJ(5)) * t30 + t58;
t63 = -qJ(3) - pkin(7);
t36 = t63 * t53;
t37 = t63 * t55;
t22 = -t51 * t36 - t49 * t37;
t13 = t32 * pkin(4) + t22;
t7 = t50 * t10 + t48 * t13;
t20 = t31 * t32;
t21 = t33 * t32;
t70 = t41 * t30;
t69 = t48 * t30;
t68 = t48 * t32;
t67 = t50 * t30;
t66 = t50 * t32;
t62 = t48 ^ 2 + t50 ^ 2;
t24 = t49 * t36 - t51 * t37;
t61 = t22 ^ 2 + t24 ^ 2;
t12 = t50 * t13;
t6 = -t48 * t10 + t12;
t3 = t7 * t48 + t6 * t50;
t60 = -t6 * t48 + t7 * t50;
t59 = -t32 * t40 + t70;
t57 = 0.2e1 * t22 * t32 - 0.2e1 * t24 * t30;
t29 = t71 * t50;
t28 = t71 * t48;
t26 = t62 * t40;
t19 = t30 * pkin(3) + t58;
t18 = t54 * t28 + t52 * t29;
t17 = -t52 * t28 + t54 * t29;
t15 = t33 * t30;
t14 = -t30 * pkin(4) + t24;
t8 = (-pkin(5) * t50 - pkin(4)) * t30 + t24;
t5 = pkin(8) * t67 + t7;
t4 = t32 * pkin(5) + t12 + (-pkin(8) * t30 - t10) * t48;
t2 = t52 * t4 + t54 * t5;
t1 = t54 * t4 - t52 * t5;
t9 = [1, 0, 0, t53 ^ 2, t53 * t72, 0, 0, 0, pkin(1) * t72, -0.2e1 * pkin(1) * t53, t57, t45 ^ 2 + t61, t57, -0.2e1 * t19 * t30, t19 * t75, t19 ^ 2 + t61, -0.2e1 * t14 * t67 + 0.2e1 * t6 * t32, 0.2e1 * t14 * t69 - 0.2e1 * t7 * t32, 0.2e1 * t60 * t30, t14 ^ 2 + t6 ^ 2 + t7 ^ 2, t16 ^ 2, t15 * t77, t32 * t77, -t15 * t75, t32 ^ 2, 0.2e1 * t1 * t32 - 0.2e1 * t8 * t15, 0.2e1 * t8 * t16 - 0.2e1 * t2 * t32; 0, 0, 0, 0, 0, t53, t55, 0, -t53 * pkin(7), -t55 * pkin(7) (-t30 * t49 - t32 * t51) * pkin(2) (-t22 * t51 + t24 * t49) * pkin(2), t44 * t32 - t70, t22, t24, t22 * t44 + t24 * t41, t14 * t48 - t59 * t50, t14 * t50 + t59 * t48, -t3, t14 * t41 + t3 * t40, t16 * t33, t33 * t15 - t16 * t31, t21, -t20, 0, -t35 * t15 + t17 * t32 + t8 * t31, t35 * t16 - t18 * t32 + t8 * t33; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t49 ^ 2 + t51 ^ 2) * pkin(2) ^ 2, 0, 0.2e1 * t44, t73, t44 ^ 2 + t76, t48 * t73, t50 * t73, -0.2e1 * t26, t62 * t40 ^ 2 + t76, t33 ^ 2, -0.2e1 * t33 * t31, 0, 0, 0, t31 * t74, t33 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, -t30, -t32, t19, -t68, -t66, t62 * t30, t60, 0, 0, 0, 0, 0, -t20, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, t22, t66, -t68, 0, t3, 0, 0, 0, 0, 0, t21, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t44, 0, 0, -t62, t26, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t69, 0, t14, 0, 0, 0, 0, 0, -t15, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t50, 0, t41, 0, 0, 0, 0, 0, t31, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t15, t32, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, 0, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
