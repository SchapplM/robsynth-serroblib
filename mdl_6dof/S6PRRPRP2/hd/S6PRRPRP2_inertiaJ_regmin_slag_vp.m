% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:47:42
% EndTime: 2019-05-05 03:47:43
% DurationCPUTime: 0.55s
% Computational Cost: add. (662->98), mult. (1381->183), div. (0->0), fcn. (1671->10), ass. (0->66)
t41 = sin(qJ(3));
t42 = sin(qJ(2));
t38 = sin(pkin(6));
t72 = cos(qJ(3));
t58 = t38 * t72;
t62 = cos(pkin(6));
t23 = t62 * t41 + t42 * t58;
t37 = sin(pkin(11));
t39 = cos(pkin(11));
t70 = t38 * t42;
t46 = -t41 * t70 + t62 * t72;
t9 = t37 * t23 - t39 * t46;
t77 = t9 ^ 2;
t43 = cos(qJ(5));
t76 = -0.2e1 * t43;
t26 = t37 * t41 - t39 * t72;
t75 = t26 * pkin(5);
t40 = sin(qJ(5));
t74 = t9 * t40;
t73 = t9 * t43;
t32 = t37 * pkin(3) + pkin(9);
t71 = t26 * t32;
t44 = cos(qJ(2));
t69 = t38 * t44;
t20 = t40 * t26;
t27 = t37 * t72 + t39 * t41;
t68 = t40 * t27;
t67 = t40 * t32;
t66 = t40 * t43;
t21 = t43 * t26;
t22 = t43 * t27;
t65 = t43 * t32;
t34 = -t72 * pkin(3) - pkin(2);
t14 = t26 * pkin(4) - t27 * pkin(9) + t34;
t60 = t72 * pkin(8);
t29 = t72 * qJ(4) + t60;
t57 = (-qJ(4) - pkin(8)) * t41;
t18 = t39 * t29 + t37 * t57;
t4 = t40 * t14 + t43 * t18;
t35 = t40 ^ 2;
t36 = t43 ^ 2;
t64 = t35 + t36;
t63 = t26 * qJ(6);
t61 = 0.2e1 * t72;
t33 = -t39 * pkin(3) - pkin(4);
t11 = t39 * t23 + t37 * t46;
t7 = t40 * t11 + t43 * t69;
t59 = -t7 * t26 + t9 * t68;
t56 = -t43 * t14 + t40 * t18;
t16 = t37 * t29 - t39 * t57;
t1 = t63 + t4;
t2 = t56 - t75;
t55 = t1 * t43 + t2 * t40;
t54 = t1 * t40 - t2 * t43;
t8 = t43 * t11 - t40 * t69;
t53 = t8 * t40 - t7 * t43;
t52 = t7 * t40 + t8 * t43;
t51 = t43 * pkin(5) + t40 * qJ(6);
t50 = pkin(5) * t40 - t43 * qJ(6);
t24 = t33 - t51;
t49 = t24 * t27 - t71;
t48 = t27 * t33 - t71;
t47 = t9 * t22 - t8 * t26;
t25 = t27 ^ 2;
t5 = t50 * t27 + t16;
t3 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 ^ 2 * t44 ^ 2 + t11 ^ 2 + t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t8 ^ 2 + t77; 0, 0, t69, -t70, 0, 0, 0, 0, 0, t44 * t58, -t41 * t69, -t11 * t26 + t9 * t27, t11 * t18 + t9 * t16 - t34 * t69, 0, 0, 0, 0, 0, t59, t47, t59, -t53 * t27, -t47, t8 * t1 + t7 * t2 + t9 * t5; 0, 1, 0, 0, t41 ^ 2, t41 * t61, 0, 0, 0, pkin(2) * t61, -0.2e1 * pkin(2) * t41, 0.2e1 * t16 * t27 - 0.2e1 * t18 * t26, t16 ^ 2 + t18 ^ 2 + t34 ^ 2, t36 * t25, -0.2e1 * t25 * t66, 0.2e1 * t26 * t22, -0.2e1 * t26 * t68, t26 ^ 2, 0.2e1 * t16 * t68 - 0.2e1 * t26 * t56, 0.2e1 * t16 * t22 - 0.2e1 * t4 * t26, -0.2e1 * t2 * t26 + 0.2e1 * t5 * t68, -0.2e1 * t54 * t27, 0.2e1 * t1 * t26 - 0.2e1 * t5 * t22, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t23, 0 (t11 * t37 - t39 * t9) * pkin(3), 0, 0, 0, 0, 0, -t73, t74, -t73, t52, -t74, t9 * t24 + t52 * t32; 0, 0, 0, 0, 0, 0, t41, t72, 0, -t41 * pkin(8), -t60 (-t26 * t37 - t27 * t39) * pkin(3) (-t16 * t39 + t18 * t37) * pkin(3), t40 * t22 (-t35 + t36) * t27, t20, t21, 0, -t16 * t43 + t48 * t40, t16 * t40 + t48 * t43, t49 * t40 - t5 * t43, t55, -t5 * t40 - t49 * t43, t5 * t24 + t55 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t37 ^ 2 + t39 ^ 2) * pkin(3) ^ 2, t35, 0.2e1 * t66, 0, 0, 0, t33 * t76, 0.2e1 * t33 * t40, t24 * t76, 0.2e1 * t64 * t32, -0.2e1 * t24 * t40, t64 * t32 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, t21, -t20, t21, -t64 * t27, t20, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -t7, 0, t8, -t7 * pkin(5) + t8 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t68, t26, -t56, -t4, -t56 + 0.2e1 * t75, -t51 * t27, 0.2e1 * t63 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, 0, -t67, -t65, -t67, -t50, t65, -t50 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t40, t43, 0, t40, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t22, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
