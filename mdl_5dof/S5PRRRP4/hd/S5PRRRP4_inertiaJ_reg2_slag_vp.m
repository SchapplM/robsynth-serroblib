% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:20
% EndTime: 2019-12-05 16:46:22
% DurationCPUTime: 0.60s
% Computational Cost: add. (190->53), mult. (418->80), div. (0->0), fcn. (408->6), ass. (0->48)
t36 = sin(qJ(4));
t34 = t36 ^ 2;
t39 = cos(qJ(4));
t35 = t39 ^ 2;
t44 = t34 + t35;
t37 = sin(qJ(3));
t60 = t37 * pkin(2);
t28 = pkin(7) + t60;
t47 = t44 * t28;
t45 = t44 * pkin(7);
t38 = sin(qJ(2));
t40 = cos(qJ(3));
t41 = cos(qJ(2));
t12 = t37 * t38 - t40 * t41;
t66 = t12 ^ 2;
t65 = -0.2e1 * t36;
t64 = -0.2e1 * t39;
t14 = t37 * t41 + t40 * t38;
t63 = t47 * t14;
t62 = t45 * t14;
t61 = t36 * pkin(7);
t59 = t39 * pkin(7);
t58 = t40 * pkin(2);
t29 = -pkin(3) - t58;
t57 = pkin(3) - t29;
t56 = t12 * t39;
t55 = t36 * t14;
t54 = t36 * t28;
t53 = t36 * t39;
t52 = t39 * t14;
t51 = t39 * t28;
t16 = -t39 * pkin(4) - t36 * qJ(5) - pkin(3);
t10 = t16 - t58;
t50 = -t10 - t16;
t49 = t47 * pkin(7);
t48 = t44 * t28 ^ 2;
t46 = t44 * pkin(7) ^ 2;
t19 = -t36 * pkin(4) + t39 * qJ(5);
t25 = -0.2e1 * t53;
t24 = 0.2e1 * t53;
t15 = 0.2e1 * t45;
t11 = t14 ^ 2;
t9 = t12 * t36;
t6 = 0.2e1 * t47;
t3 = t45 + t47;
t2 = t44 * t14;
t1 = t44 * t11 + t66;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 ^ 2 + t41 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 + t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t38, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t14, 0, (-t12 * t40 + t14 * t37) * pkin(2), 0, 0, 0, 0, 0, 0, -t56, t9, t2, t12 * t29 + t63, 0, 0, 0, 0, 0, 0, -t56, t2, -t9, t12 * t10 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t60, 0, (t37 ^ 2 + t40 ^ 2) * pkin(2) ^ 2, t34, t24, 0, t35, 0, 0, t29 * t64, 0.2e1 * t29 * t36, t6, t29 ^ 2 + t48, t34, 0, t25, 0, 0, t35, t10 * t64, t6, t10 * t65, t10 ^ 2 + t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t14, 0, 0, 0, 0, 0, 0, 0, 0, -t56, t9, t2, -t12 * pkin(3) + t62, 0, 0, 0, 0, 0, 0, -t56, t2, -t9, t12 * t16 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t58, -t60, 0, 0, t34, t24, 0, t35, 0, 0, t57 * t39, -t57 * t36, t3, -t29 * pkin(3) + t49, t34, 0, t25, 0, 0, t35, t50 * t39, t3, t50 * t36, t10 * t16 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t34, t24, 0, t35, 0, 0, 0.2e1 * pkin(3) * t39, pkin(3) * t65, t15, pkin(3) ^ 2 + t46, t34, 0, t25, 0, 0, t35, t16 * t64, t15, t16 * t65, t16 ^ 2 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t52, 0, 0, 0, 0, 0, 0, 0, 0, -t55, 0, t52, t19 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t39, 0, -t54, -t51, 0, 0, 0, t36, 0, 0, -t39, 0, -t54, t19, t51, t19 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t39, 0, -t61, -t59, 0, 0, 0, t36, 0, 0, -t39, 0, -t61, t19, t59, t19 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
