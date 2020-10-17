% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPP6
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:00
% EndTime: 2019-12-31 21:02:03
% DurationCPUTime: 0.92s
% Computational Cost: add. (709->116), mult. (1499->227), div. (0->0), fcn. (1571->6), ass. (0->77)
t48 = sin(pkin(8));
t49 = cos(pkin(8));
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t29 = t48 * t52 + t49 * t50;
t51 = sin(qJ(2));
t23 = t29 * t51;
t93 = t23 ^ 2;
t27 = t48 * t50 - t49 * t52;
t92 = t27 ^ 2;
t41 = -t52 * pkin(3) - pkin(2);
t91 = 0.2e1 * t41;
t90 = -0.2e1 * t51;
t53 = cos(qJ(2));
t89 = 0.2e1 * t53;
t88 = pkin(2) * t52;
t87 = pkin(6) * t50;
t45 = t51 ^ 2;
t86 = t45 * pkin(6);
t85 = t48 * pkin(3);
t84 = t49 * pkin(3);
t43 = t51 * pkin(6);
t83 = t53 * pkin(3);
t33 = -t53 * pkin(2) - t51 * pkin(7) - pkin(1);
t30 = t52 * t33;
t67 = qJ(4) * t51;
t12 = -t52 * t67 + t30 + (-pkin(3) - t87) * t53;
t72 = t52 * t53;
t64 = pkin(6) * t72;
t14 = t64 + (t33 - t67) * t50;
t4 = t48 * t12 + t49 * t14;
t69 = -qJ(4) - pkin(7);
t34 = t69 * t52;
t61 = t69 * t50;
t16 = -t48 * t34 - t49 * t61;
t82 = t16 * t53;
t18 = -t49 * t34 + t48 * t61;
t81 = t18 * t53;
t80 = t23 * t27;
t73 = t52 * t51;
t76 = t50 * t51;
t25 = -t48 * t76 + t49 * t73;
t79 = t25 * t23;
t78 = t29 * t27;
t77 = t29 * t53;
t75 = t50 * t52;
t74 = t50 * t53;
t71 = t53 * t23;
t70 = t53 * t27;
t32 = pkin(3) * t76 + t43;
t44 = t50 ^ 2;
t46 = t52 ^ 2;
t68 = t44 + t46;
t66 = t51 * t89;
t65 = t16 ^ 2 + t18 ^ 2;
t63 = t50 * t73;
t62 = t16 * t25 - t18 * t23;
t60 = -t49 * t12 + t48 * t14;
t20 = -pkin(6) * t74 + t30;
t21 = t50 * t33 + t64;
t59 = -t20 * t50 + t21 * t52;
t58 = -t29 * t23 - t25 * t27;
t57 = 0.2e1 * t16 * t29 - 0.2e1 * t18 * t27;
t55 = pkin(6) ^ 2;
t47 = t53 ^ 2;
t42 = t45 * t55;
t39 = pkin(4) + t84;
t37 = qJ(5) + t85;
t26 = t29 ^ 2;
t22 = t25 ^ 2;
t19 = -0.2e1 * t25 * t53;
t11 = t25 * t29;
t8 = t27 * pkin(4) - t29 * qJ(5) + t41;
t5 = t23 * pkin(4) - t25 * qJ(5) + t32;
t2 = t53 * pkin(4) + t60;
t1 = -t53 * qJ(5) + t4;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t45, t66, 0, t47, 0, 0, pkin(1) * t89, pkin(1) * t90, 0.2e1 * (t45 + t47) * pkin(6), pkin(1) ^ 2 + t47 * t55 + t42, t46 * t45, -0.2e1 * t45 * t75, t72 * t90, t44 * t45, t50 * t66, t47, -0.2e1 * t20 * t53 + 0.2e1 * t50 * t86, 0.2e1 * t21 * t53 + 0.2e1 * t52 * t86, 0.2e1 * (-t20 * t52 - t21 * t50) * t51, t20 ^ 2 + t21 ^ 2 + t42, t22, -0.2e1 * t79, t19, t93, 0.2e1 * t71, t47, 0.2e1 * t32 * t23 + 0.2e1 * t53 * t60, 0.2e1 * t32 * t25 + 0.2e1 * t4 * t53, -0.2e1 * t4 * t23 + 0.2e1 * t25 * t60, t32 ^ 2 + t4 ^ 2 + t60 ^ 2, t22, t19, 0.2e1 * t79, t47, -0.2e1 * t71, t93, 0.2e1 * t2 * t53 + 0.2e1 * t5 * t23, -0.2e1 * t1 * t23 + 0.2e1 * t2 * t25, -0.2e1 * t1 * t53 - 0.2e1 * t5 * t25, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t53, 0, -t43, -t53 * pkin(6), 0, 0, t63, (-t44 + t46) * t51, -t74, -t63, -t72, 0, -pkin(6) * t73 + (-pkin(2) * t51 + pkin(7) * t53) * t50, pkin(7) * t72 + (t87 - t88) * t51, t59, -pkin(2) * t43 + pkin(7) * t59, t11, t58, -t77, t80, t70, 0, t41 * t23 + t32 * t27 + t82, t41 * t25 + t32 * t29 + t81, -t4 * t27 + t29 * t60 + t62, t16 * t60 + t4 * t18 + t32 * t41, t11, -t77, -t58, 0, -t70, t80, t8 * t23 + t5 * t27 + t82, -t1 * t27 + t2 * t29 + t62, -t8 * t25 - t5 * t29 - t81, t1 * t18 + t2 * t16 + t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t44, 0.2e1 * t75, 0, t46, 0, 0, 0.2e1 * t88, -0.2e1 * pkin(2) * t50, 0.2e1 * t68 * pkin(7), pkin(7) ^ 2 * t68 + pkin(2) ^ 2, t26, -0.2e1 * t78, 0, t92, 0, 0, t27 * t91, t29 * t91, t57, t41 ^ 2 + t65, t26, 0, 0.2e1 * t78, 0, 0, t92, 0.2e1 * t8 * t27, t57, -0.2e1 * t8 * t29, t8 ^ 2 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, -t76, -t53, t20, -t21, 0, 0, 0, 0, t25, 0, -t23, -t53, -t49 * t83 - t60, t48 * t83 - t4, (-t23 * t48 - t25 * t49) * pkin(3), (t4 * t48 - t49 * t60) * pkin(3), 0, t25, 0, -t53, t23, 0, (-pkin(4) - t39) * t53 - t60, -t37 * t23 - t39 * t25, (-qJ(5) - t37) * t53 + t4, t1 * t37 - t2 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t52, 0, -t50 * pkin(7), -t52 * pkin(7), 0, 0, 0, 0, t29, 0, -t27, 0, -t16, -t18, (-t27 * t48 - t29 * t49) * pkin(3), (-t16 * t49 + t18 * t48) * pkin(3), 0, t29, 0, 0, t27, 0, -t16, -t37 * t27 - t39 * t29, t18, -t16 * t39 + t18 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t84, -0.2e1 * t85, 0, (t48 ^ 2 + t49 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t39, 0, 0.2e1 * t37, t37 ^ 2 + t39 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t25, 0, t32, 0, 0, 0, 0, 0, 0, t23, 0, -t25, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t29, 0, t41, 0, 0, 0, 0, 0, 0, t27, 0, -t29, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t25, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
