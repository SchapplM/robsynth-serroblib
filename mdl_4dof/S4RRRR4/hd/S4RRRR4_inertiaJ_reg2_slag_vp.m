% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:15
% EndTime: 2019-12-31 17:26:16
% DurationCPUTime: 0.41s
% Computational Cost: add. (331->59), mult. (701->131), div. (0->0), fcn. (730->6), ass. (0->50)
t34 = sin(qJ(3));
t56 = t34 * pkin(2);
t24 = pkin(7) + t56;
t33 = sin(qJ(4));
t29 = t33 ^ 2;
t36 = cos(qJ(4));
t31 = t36 ^ 2;
t47 = t29 + t31;
t49 = t47 * t24;
t38 = cos(qJ(2));
t57 = -pkin(6) - pkin(5);
t18 = t57 * t38;
t37 = cos(qJ(3));
t35 = sin(qJ(2));
t44 = t57 * t35;
t7 = -t34 * t18 - t37 * t44;
t61 = t7 ^ 2;
t14 = t34 * t35 - t37 * t38;
t60 = t14 ^ 2;
t26 = -t38 * pkin(2) - pkin(1);
t59 = 0.2e1 * t26;
t58 = 0.2e1 * t38;
t55 = t37 * pkin(2);
t54 = t7 * t36;
t25 = -pkin(3) - t55;
t53 = pkin(3) - t25;
t16 = t34 * t38 + t37 * t35;
t52 = t33 * t16;
t51 = t33 * t36;
t50 = t36 * t16;
t48 = t47 * pkin(7);
t30 = t35 ^ 2;
t32 = t38 ^ 2;
t46 = t30 + t32;
t45 = -0.2e1 * t16 * t14;
t43 = -pkin(3) * t16 - pkin(7) * t14;
t4 = t14 * pkin(3) - t16 * pkin(7) + t26;
t9 = -t37 * t18 + t34 * t44;
t2 = -t33 * t9 + t36 * t4;
t3 = t33 * t4 + t36 * t9;
t1 = -t2 * t33 + t3 * t36;
t42 = -t14 * t24 + t16 * t25;
t21 = 0.2e1 * t51;
t13 = t16 ^ 2;
t12 = t36 * t14;
t11 = t33 * t14;
t10 = t33 * t50;
t6 = t7 * t33;
t5 = (-t29 + t31) * t16;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t30, t35 * t58, 0, t32, 0, 0, pkin(1) * t58, -0.2e1 * pkin(1) * t35, 0.2e1 * t46 * pkin(5), t46 * pkin(5) ^ 2 + pkin(1) ^ 2, t13, t45, 0, t60, 0, 0, t14 * t59, t16 * t59, -0.2e1 * t9 * t14 + 0.2e1 * t7 * t16, t26 ^ 2 + t9 ^ 2 + t61, t31 * t13, -0.2e1 * t13 * t51, 0.2e1 * t14 * t50, t29 * t13, t33 * t45, t60, 0.2e1 * t2 * t14 + 0.2e1 * t7 * t52, -0.2e1 * t3 * t14 + 0.2e1 * t7 * t50, 0.2e1 * (-t2 * t36 - t3 * t33) * t16, t2 ^ 2 + t3 ^ 2 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t38, 0, -t35 * pkin(5), -t38 * pkin(5), 0, 0, 0, 0, t16, 0, -t14, 0, -t7, -t9, (-t14 * t34 - t16 * t37) * pkin(2), (t34 * t9 - t37 * t7) * pkin(2), t10, t5, t11, -t10, t12, 0, t42 * t33 - t54, t42 * t36 + t6, t1, t1 * t24 + t7 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t55, -0.2e1 * t56, 0, (t34 ^ 2 + t37 ^ 2) * pkin(2) ^ 2, t29, t21, 0, t31, 0, 0, -0.2e1 * t25 * t36, 0.2e1 * t25 * t33, 0.2e1 * t49, t47 * t24 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, -t7, -t9, 0, 0, t10, t5, t11, -t10, t12, 0, t43 * t33 - t54, t43 * t36 + t6, t1, -t7 * pkin(3) + t1 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t55, -t56, 0, 0, t29, t21, 0, t31, 0, 0, t53 * t36, -t53 * t33, t48 + t49, -t25 * pkin(3) + pkin(7) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t29, t21, 0, t31, 0, 0, 0.2e1 * pkin(3) * t36, -0.2e1 * pkin(3) * t33, 0.2e1 * t48, t47 * pkin(7) ^ 2 + pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t52, t14, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, t36, 0, -t33 * t24, -t36 * t24, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, t36, 0, -t33 * pkin(7), -t36 * pkin(7), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
