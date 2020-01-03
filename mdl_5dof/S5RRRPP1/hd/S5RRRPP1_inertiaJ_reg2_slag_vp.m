% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t52 = sin(qJ(2));
t75 = t52 * pkin(1);
t41 = pkin(7) + t75;
t51 = sin(qJ(3));
t47 = t51 ^ 2;
t53 = cos(qJ(3));
t48 = t53 ^ 2;
t63 = t47 + t48;
t65 = t63 * t41;
t44 = t53 * qJ(4);
t74 = t53 * pkin(7);
t32 = t44 + t74;
t49 = sin(pkin(8));
t50 = cos(pkin(8));
t59 = (-qJ(4) - pkin(7)) * t51;
t19 = t49 * t32 - t50 * t59;
t21 = t50 * t32 + t49 * t59;
t27 = t49 * t51 - t50 * t53;
t29 = t49 * t53 + t50 * t51;
t88 = t19 * t29 - t21 * t27;
t92 = 0.2e1 * t88;
t67 = t53 * t41;
t23 = t44 + t67;
t58 = (-qJ(4) - t41) * t51;
t11 = t49 * t23 - t50 * t58;
t13 = t50 * t23 + t49 * t58;
t89 = t11 * t29 - t13 * t27;
t91 = 0.2e1 * t89;
t90 = t88 + t89;
t87 = t27 ^ 2;
t86 = 0.2e1 * t27;
t85 = -0.2e1 * t29;
t43 = -pkin(3) * t53 - pkin(2);
t54 = cos(qJ(2));
t73 = t54 * pkin(1);
t31 = t43 - t73;
t84 = 0.2e1 * t31;
t83 = 0.2e1 * t43;
t82 = 0.2e1 * t53;
t77 = t49 * pkin(3);
t76 = t50 * pkin(3);
t42 = -pkin(2) - t73;
t72 = pkin(2) - t42;
t15 = pkin(4) * t27 - qJ(5) * t29 + t43;
t9 = t15 - t73;
t71 = t15 + t9;
t68 = t29 * t27;
t66 = t31 + t43;
t64 = t63 * pkin(7);
t62 = t11 ^ 2 + t13 ^ 2;
t61 = t19 ^ 2 + t21 ^ 2;
t60 = t11 * t19 + t13 * t21;
t39 = pkin(4) + t76;
t36 = qJ(5) + t77;
t35 = t51 * t82;
t26 = t29 ^ 2;
t17 = -0.2e1 * t68;
t16 = 0.2e1 * t68;
t14 = (-t27 * t49 - t29 * t50) * pkin(3);
t5 = -t27 * t36 - t29 * t39;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t73, -0.2e1 * t75, 0, (t52 ^ 2 + t54 ^ 2) * pkin(1) ^ 2, t47, t35, 0, t48, 0, 0, -0.2e1 * t42 * t53, 0.2e1 * t42 * t51, 0.2e1 * t65, t41 ^ 2 * t63 + t42 ^ 2, t26, t17, 0, t87, 0, 0, t27 * t84, t29 * t84, t91, t31 ^ 2 + t62, t26, 0, t16, 0, 0, t87, t9 * t86, t91, t9 * t85, t9 ^ 2 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t73, -t75, 0, 0, t47, t35, 0, t48, 0, 0, t72 * t53, -t72 * t51, t64 + t65, -t42 * pkin(2) + pkin(7) * t65, t26, t17, 0, t87, 0, 0, t66 * t27, t66 * t29, t90, t31 * t43 + t60, t26, 0, t16, 0, 0, t87, t71 * t27, t90, -t71 * t29, t15 * t9 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t47, t35, 0, t48, 0, 0, pkin(2) * t82, -0.2e1 * pkin(2) * t51, 0.2e1 * t64, pkin(7) ^ 2 * t63 + pkin(2) ^ 2, t26, t17, 0, t87, 0, 0, t27 * t83, t29 * t83, t92, t43 ^ 2 + t61, t26, 0, t16, 0, 0, t87, t15 * t86, t92, t15 * t85, t15 ^ 2 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t53, 0, -t51 * t41, -t67, 0, 0, 0, 0, t29, 0, -t27, 0, -t11, -t13, t14, (-t11 * t50 + t13 * t49) * pkin(3), 0, t29, 0, 0, t27, 0, -t11, t5, t13, -t11 * t39 + t13 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t53, 0, -t51 * pkin(7), -t74, 0, 0, 0, 0, t29, 0, -t27, 0, -t19, -t21, t14, (-t19 * t50 + t21 * t49) * pkin(3), 0, t29, 0, 0, t27, 0, -t19, t5, t21, -t19 * t39 + t21 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t76, -0.2e1 * t77, 0, (t49 ^ 2 + t50 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t39, 0, 0.2e1 * t36, t36 ^ 2 + t39 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t29, 0, t31, 0, 0, 0, 0, 0, 0, t27, 0, -t29, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t29, 0, t43, 0, 0, 0, 0, 0, 0, t27, 0, -t29, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
