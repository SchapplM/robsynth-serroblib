% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t53 = sin(qJ(2));
t80 = t53 * pkin(1);
t40 = pkin(7) + t80;
t52 = sin(qJ(3));
t49 = t52 ^ 2;
t55 = cos(qJ(3));
t50 = t55 ^ 2;
t67 = t49 + t50;
t69 = t67 * t40;
t48 = t55 * pkin(8);
t78 = t55 * pkin(7);
t33 = t48 + t78;
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t64 = (-pkin(8) - pkin(7)) * t52;
t20 = t51 * t33 - t54 * t64;
t22 = t54 * t33 + t51 * t64;
t28 = t51 * t52 - t54 * t55;
t30 = t51 * t55 + t54 * t52;
t91 = t20 * t30 - t22 * t28;
t95 = 0.2e1 * t91;
t71 = t55 * t40;
t24 = t48 + t71;
t63 = (-pkin(8) - t40) * t52;
t11 = t51 * t24 - t54 * t63;
t13 = t54 * t24 + t51 * t63;
t92 = t11 * t30 - t13 * t28;
t94 = 0.2e1 * t92;
t93 = t91 + t92;
t90 = t28 ^ 2;
t89 = 0.2e1 * t28;
t88 = -0.2e1 * t30;
t44 = -t55 * pkin(3) - pkin(2);
t56 = cos(qJ(2));
t77 = t56 * pkin(1);
t32 = t44 - t77;
t87 = 0.2e1 * t32;
t86 = 0.2e1 * t44;
t85 = 0.2e1 * t55;
t79 = t54 * pkin(3);
t43 = -pkin(2) - t77;
t76 = pkin(2) - t43;
t14 = t28 * pkin(4) - t30 * qJ(5) + t44;
t9 = t14 - t77;
t75 = t14 + t9;
t72 = t30 * t28;
t70 = t32 + t44;
t68 = t67 * pkin(7);
t66 = t11 ^ 2 + t13 ^ 2;
t65 = t20 ^ 2 + t22 ^ 2;
t62 = t11 * t20 + t13 * t22;
t58 = 0.2e1 * pkin(4);
t57 = 0.2e1 * qJ(5);
t45 = t51 * pkin(3);
t41 = pkin(4) + t79;
t38 = t45 + qJ(5);
t36 = t52 * t85;
t27 = t30 ^ 2;
t18 = -0.2e1 * t72;
t17 = 0.2e1 * t72;
t16 = -pkin(4) * t30 - t28 * qJ(5);
t15 = (-t28 * t51 - t30 * t54) * pkin(3);
t5 = -t38 * t28 - t41 * t30;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t77, -0.2e1 * t80, 0, (t53 ^ 2 + t56 ^ 2) * pkin(1) ^ 2, t49, t36, 0, t50, 0, 0, -0.2e1 * t43 * t55, 0.2e1 * t43 * t52, 0.2e1 * t69, t67 * t40 ^ 2 + t43 ^ 2, t27, t18, 0, t90, 0, 0, t28 * t87, t30 * t87, t94, t32 ^ 2 + t66, t27, 0, t17, 0, 0, t90, t9 * t89, t94, t9 * t88, t9 ^ 2 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t77, -t80, 0, 0, t49, t36, 0, t50, 0, 0, t76 * t55, -t76 * t52, t68 + t69, -t43 * pkin(2) + pkin(7) * t69, t27, t18, 0, t90, 0, 0, t70 * t28, t70 * t30, t93, t32 * t44 + t62, t27, 0, t17, 0, 0, t90, t75 * t28, t93, -t75 * t30, t9 * t14 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t49, t36, 0, t50, 0, 0, pkin(2) * t85, -0.2e1 * pkin(2) * t52, 0.2e1 * t68, t67 * pkin(7) ^ 2 + pkin(2) ^ 2, t27, t18, 0, t90, 0, 0, t28 * t86, t30 * t86, t95, t44 ^ 2 + t65, t27, 0, t17, 0, 0, t90, t14 * t89, t95, t14 * t88, t14 ^ 2 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t55, 0, -t52 * t40, -t71, 0, 0, 0, 0, t30, 0, -t28, 0, -t11, -t13, t15, (-t11 * t54 + t13 * t51) * pkin(3), 0, t30, 0, 0, t28, 0, -t11, t5, t13, -t11 * t41 + t13 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t55, 0, -t52 * pkin(7), -t78, 0, 0, 0, 0, t30, 0, -t28, 0, -t20, -t22, t15, (-t20 * t54 + t22 * t51) * pkin(3), 0, t30, 0, 0, t28, 0, -t20, t5, t22, -t20 * t41 + t22 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t79, -0.2e1 * t45, 0, (t51 ^ 2 + t54 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t41, 0, 0.2e1 * t38, t38 ^ 2 + t41 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, -t28, 0, -t11, -t13, 0, 0, 0, t30, 0, 0, t28, 0, -t11, t16, t13, -t11 * pkin(4) + t13 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, -t28, 0, -t20, -t22, 0, 0, 0, t30, 0, 0, t28, 0, -t20, t16, t22, -t20 * pkin(4) + t22 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t79, -t45, 0, 0, 0, 0, 0, 1, 0, 0, t58 + t79, 0, t57 + t45, t41 * pkin(4) + t38 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t58, 0, t57, pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
