% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t53 = sin(pkin(10));
t54 = cos(pkin(10));
t57 = sin(qJ(2));
t73 = cos(qJ(2));
t34 = t53 * t57 - t54 * t73;
t80 = -0.2e1 * t34;
t79 = 0.2e1 * t34;
t46 = -t54 * pkin(2) - pkin(3);
t58 = cos(qJ(4));
t41 = -t58 * pkin(4) + t46;
t78 = 0.2e1 * t41;
t77 = t34 * pkin(4);
t55 = sin(qJ(5));
t56 = sin(qJ(4));
t72 = cos(qJ(5));
t61 = t72 * t58;
t38 = t55 * t56 - t61;
t76 = t38 * pkin(5);
t75 = t55 * pkin(4);
t45 = t53 * pkin(2) + pkin(8);
t74 = pkin(9) + t45;
t35 = t53 * t73 + t54 * t57;
t68 = t56 * t35;
t16 = t35 * t61 - t55 * t68;
t71 = t16 * t38;
t40 = t55 * t58 + t56 * t72;
t70 = t40 * t34;
t69 = t56 * t34;
t67 = t56 * t58;
t42 = (-qJ(3) - pkin(7)) * t57;
t63 = t73 * pkin(7);
t43 = qJ(3) * t73 + t63;
t23 = t53 * t42 + t54 * t43;
t66 = t58 * t23;
t65 = t58 * t35;
t64 = 0.2e1 * t73;
t50 = t72 * pkin(4);
t49 = -pkin(2) * t73 - pkin(1);
t20 = t34 * pkin(3) - t35 * pkin(8) + t49;
t7 = t66 + (-pkin(9) * t35 + t20) * t56;
t62 = t72 * t7;
t9 = t58 * t20 - t56 * t23;
t6 = -pkin(9) * t65 + t77 + t9;
t3 = -t55 * t7 + t72 * t6;
t30 = t74 * t56;
t31 = t74 * t58;
t18 = -t72 * t30 - t55 * t31;
t21 = -t54 * t42 + t53 * t43;
t14 = pkin(4) * t68 + t21;
t60 = -t34 * t45 + t35 * t46;
t4 = t55 * t6 + t62;
t19 = -t55 * t30 + t31 * t72;
t52 = t58 ^ 2;
t51 = t56 ^ 2;
t48 = t50 + pkin(5);
t36 = t40 ^ 2;
t33 = t35 ^ 2;
t32 = t34 ^ 2;
t29 = t58 * t34;
t25 = t41 + t76;
t24 = t38 * t34;
t15 = t40 * t35;
t13 = -t38 * qJ(6) + t19;
t12 = -t40 * qJ(6) + t18;
t11 = t40 * t15;
t10 = t56 * t20 + t66;
t8 = t15 * pkin(5) + t14;
t2 = -t15 * qJ(6) + t4;
t1 = t34 * pkin(5) - t16 * qJ(6) + t3;
t5 = [1, 0, 0, t57 ^ 2, t57 * t64, 0, 0, 0, pkin(1) * t64, -0.2e1 * pkin(1) * t57, 0.2e1 * t21 * t35 - 0.2e1 * t23 * t34, t21 ^ 2 + t23 ^ 2 + t49 ^ 2, t52 * t33, -0.2e1 * t33 * t67, t65 * t79, t68 * t80, t32, 0.2e1 * t21 * t68 + 0.2e1 * t9 * t34, -0.2e1 * t10 * t34 + 0.2e1 * t21 * t65, t16 ^ 2, -0.2e1 * t16 * t15, t16 * t79, t15 * t80, t32, 0.2e1 * t14 * t15 + 0.2e1 * t3 * t34, 0.2e1 * t14 * t16 - 0.2e1 * t4 * t34, -0.2e1 * t1 * t16 - 0.2e1 * t2 * t15, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, t57, t73, 0, -t57 * pkin(7), -t63 (-t34 * t53 - t35 * t54) * pkin(2) (-t21 * t54 + t23 * t53) * pkin(2), t56 * t65 (-t51 + t52) * t35, t69, t29, 0, -t21 * t58 + t56 * t60, t21 * t56 + t58 * t60, t16 * t40, -t11 - t71, t70, -t24, 0, t14 * t38 + t41 * t15 + t18 * t34, t14 * t40 + t41 * t16 - t19 * t34, -t1 * t40 - t12 * t16 - t13 * t15 - t2 * t38, t1 * t12 + t2 * t13 + t8 * t25; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t53 ^ 2 + t54 ^ 2) * pkin(2) ^ 2, t51, 0.2e1 * t67, 0, 0, 0, -0.2e1 * t46 * t58, 0.2e1 * t46 * t56, t36, -0.2e1 * t40 * t38, 0, 0, 0, t38 * t78, t40 * t78, -0.2e1 * t12 * t40 - 0.2e1 * t13 * t38, t12 ^ 2 + t13 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, t29, -t69, 0, 0, 0, 0, 0, -t24, -t70, -t11 + t71, -t1 * t38 + t2 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t38 + t13 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 ^ 2 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t68, t34, t9, -t10, 0, 0, t16, -t15, t34, t34 * t50 + t3, -t62 + (-t6 - t77) * t55, -t15 * t75 - t48 * t16, t1 * t48 + t2 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t58, 0, -t56 * t45, -t58 * t45, 0, 0, t40, -t38, 0, t18, -t19, -t38 * t75 - t48 * t40, t12 * t48 + t13 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t56, 0, 0, 0, 0, 0, -t38, -t40, 0, -t38 * t48 + t40 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t50, -0.2e1 * t75, 0, t55 ^ 2 * pkin(4) ^ 2 + t48 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t34, t3, -t4, -pkin(5) * t16, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t38, 0, t18, -t19, -pkin(5) * t40, t12 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t40, 0, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t50, -t75, 0, t48 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t5;
