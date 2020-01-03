% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t49 = sin(qJ(4));
t47 = t49 ^ 2;
t52 = cos(qJ(4));
t48 = t52 ^ 2;
t86 = t47 + t48;
t54 = cos(qJ(2));
t45 = t54 * pkin(1);
t39 = t45 + pkin(2);
t50 = sin(qJ(3));
t53 = cos(qJ(3));
t51 = sin(qJ(2));
t78 = t51 * pkin(1);
t58 = t53 * t78;
t20 = t39 * t50 + t58;
t18 = pkin(8) + t20;
t66 = t86 * t18;
t79 = t50 * pkin(2);
t37 = pkin(8) + t79;
t62 = t86 * t37;
t85 = -0.2e1 * t49;
t84 = -0.2e1 * t52;
t83 = t62 * t18;
t82 = t66 * pkin(8);
t81 = pkin(3) * t49;
t80 = t49 * pkin(8);
t77 = t52 * pkin(8);
t23 = -pkin(4) * t52 - qJ(5) * t49 - pkin(3);
t44 = t53 * pkin(2);
t21 = t23 - t44;
t61 = -t39 * t53 + t50 * t78;
t7 = t23 + t61;
t76 = -t21 - t7;
t75 = -t23 - t7;
t17 = -pkin(3) + t61;
t74 = t17 * t52;
t38 = -t44 - pkin(3);
t73 = t38 * t52;
t72 = t49 * t18;
t71 = t49 * t37;
t70 = t49 * t52;
t69 = t52 * t18;
t68 = t52 * t37;
t67 = t86 * t18 ^ 2;
t65 = -t21 - t23;
t64 = t62 * pkin(8);
t63 = t86 * t37 ^ 2;
t60 = t86 * pkin(8) ^ 2;
t59 = t86 * pkin(8);
t26 = -pkin(4) * t49 + qJ(5) * t52;
t46 = pkin(3) * t52;
t35 = -0.2e1 * t70;
t34 = 0.2e1 * t70;
t32 = t38 * t49;
t22 = 0.2e1 * t59;
t15 = t17 * t49;
t10 = 0.2e1 * t62;
t4 = t59 + t62;
t3 = 0.2e1 * t66;
t2 = t59 + t66;
t1 = t62 + t66;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t45, -0.2e1 * t78, 0, (t51 ^ 2 + t54 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t61, -0.2e1 * t20, 0, t20 ^ 2 + t61 ^ 2, t47, t34, 0, t48, 0, 0, -0.2e1 * t74, 0.2e1 * t15, t3, t17 ^ 2 + t67, t47, 0, t35, 0, 0, t48, t7 * t84, t3, t7 * t85, t7 ^ 2 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t45, -t78, 0, 0, 0, 0, 0, 0, 0, 1, t44 - t61, -t58 + (-pkin(2) - t39) * t50, 0, (t20 * t50 - t53 * t61) * pkin(2), t47, t34, 0, t48, 0, 0, (-t17 - t38) * t52, t32 + t15, t1, t17 * t38 + t83, t47, 0, t35, 0, 0, t48, t76 * t52, t1, t76 * t49, t21 * t7 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t44, -0.2e1 * t79, 0, (t50 ^ 2 + t53 ^ 2) * pkin(2) ^ 2, t47, t34, 0, t48, 0, 0, -0.2e1 * t73, 0.2e1 * t32, t10, t38 ^ 2 + t63, t47, 0, t35, 0, 0, t48, t21 * t84, t10, t21 * t85, t21 ^ 2 + t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t61, -t20, 0, 0, t47, t34, 0, t48, 0, 0, t46 - t74, t15 - t81, t2, -pkin(3) * t17 + t82, t47, 0, t35, 0, 0, t48, t75 * t52, t2, t75 * t49, t23 * t7 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t44, -t79, 0, 0, t47, t34, 0, t48, 0, 0, t46 - t73, t32 - t81, t4, -pkin(3) * t38 + t64, t47, 0, t35, 0, 0, t48, t65 * t52, t4, t65 * t49, t21 * t23 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t47, t34, 0, t48, 0, 0, 0.2e1 * t46, -0.2e1 * t81, t22, pkin(3) ^ 2 + t60, t47, 0, t35, 0, 0, t48, t23 * t84, t22, t23 * t85, t23 ^ 2 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t52, 0, -t72, -t69, 0, 0, 0, t49, 0, 0, -t52, 0, -t72, t26, t69, t26 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t52, 0, -t71, -t68, 0, 0, 0, t49, 0, 0, -t52, 0, -t71, t26, t68, t26 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t52, 0, -t80, -t77, 0, 0, 0, t49, 0, 0, -t52, 0, -t80, t26, t77, t26 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t5;
