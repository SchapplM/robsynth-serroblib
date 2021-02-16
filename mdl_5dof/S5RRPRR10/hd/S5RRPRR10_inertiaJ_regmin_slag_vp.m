% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:02:43
% EndTime: 2021-01-15 22:02:47
% DurationCPUTime: 0.60s
% Computational Cost: add. (687->97), mult. (1814->215), div. (0->0), fcn. (2062->10), ass. (0->78)
t42 = sin(pkin(10));
t43 = sin(pkin(5));
t44 = cos(pkin(10));
t48 = sin(qJ(2));
t51 = cos(qJ(2));
t28 = (t42 * t51 + t44 * t48) * t43;
t45 = cos(pkin(5));
t47 = sin(qJ(4));
t50 = cos(qJ(4));
t19 = t47 * t28 - t45 * t50;
t86 = -0.2e1 * t19;
t85 = 0.2e1 * t47;
t84 = pkin(1) * t48;
t49 = cos(qJ(5));
t83 = pkin(4) * t49;
t72 = t43 * t51;
t73 = t43 * t48;
t27 = t42 * t73 - t44 * t72;
t32 = (-pkin(2) * t51 - pkin(1)) * t43;
t14 = t27 * pkin(3) - t28 * pkin(8) + t32;
t34 = t45 * t51 * pkin(1);
t59 = pkin(7) + qJ(3);
t78 = t45 * pkin(2);
t21 = -t59 * t73 + t34 + t78;
t55 = t45 * t84;
t24 = t59 * t72 + t55;
t11 = t42 * t21 + t44 * t24;
t9 = t45 * pkin(8) + t11;
t6 = t50 * t14 - t47 * t9;
t3 = -t27 * pkin(4) - t6;
t46 = sin(qJ(5));
t82 = t3 * t46;
t81 = t3 * t49;
t80 = t42 * pkin(2);
t79 = t44 * pkin(2);
t20 = t50 * t28 + t45 * t47;
t13 = t49 * t20 + t46 * t27;
t77 = t13 * t46;
t76 = t13 * t50;
t36 = pkin(8) + t80;
t75 = t36 * t46;
t38 = t43 ^ 2;
t74 = t38 * t51;
t71 = t46 * t19;
t70 = t46 * t47;
t69 = t46 * t49;
t68 = t46 * t50;
t67 = t47 * t19;
t66 = t47 * t27;
t65 = t47 * t36;
t64 = t49 * t19;
t63 = t49 * t47;
t62 = t49 * t50;
t12 = t46 * t20 - t49 * t27;
t61 = t50 * t12;
t60 = t50 * t36;
t58 = -t44 * t21 + t42 * t24;
t57 = 0.2e1 * t43 * t45;
t56 = t50 * t85;
t54 = t46 * t67;
t53 = t19 * t63;
t37 = -pkin(3) - t79;
t8 = -t45 * pkin(3) + t58;
t7 = t47 * t14 + t50 * t9;
t41 = t49 ^ 2;
t40 = t47 ^ 2;
t39 = t46 ^ 2;
t31 = -t50 * pkin(4) - t47 * pkin(9) + t37;
t30 = pkin(7) * t72 + t55;
t29 = -pkin(7) * t73 + t34;
t26 = t50 * t27;
t18 = t46 * t31 + t49 * t60;
t17 = t49 * t31 - t46 * t60;
t5 = t19 * pkin(4) - t20 * pkin(9) + t8;
t4 = t27 * pkin(9) + t7;
t2 = t49 * t4 + t46 * t5;
t1 = -t46 * t4 + t49 * t5;
t10 = [1, 0, 0, t38 * t48 ^ 2, 0.2e1 * t48 * t74, t48 * t57, t51 * t57, t45 ^ 2, 0.2e1 * pkin(1) * t74 + 0.2e1 * t29 * t45, -0.2e1 * t30 * t45 - 0.2e1 * t38 * t84, 0.2e1 * t32 * t27 - 0.2e1 * t45 * t58, -0.2e1 * t11 * t45 + 0.2e1 * t32 * t28, -0.2e1 * t11 * t27 + 0.2e1 * t28 * t58, t11 ^ 2 + t32 ^ 2 + t58 ^ 2, t20 ^ 2, t20 * t86, 0.2e1 * t20 * t27, t27 * t86, t27 ^ 2, 0.2e1 * t8 * t19 + 0.2e1 * t6 * t27, 0.2e1 * t8 * t20 - 0.2e1 * t7 * t27, t13 ^ 2, -0.2e1 * t13 * t12, 0.2e1 * t13 * t19, t12 * t86, t19 ^ 2, 0.2e1 * t1 * t19 + 0.2e1 * t3 * t12, 0.2e1 * t3 * t13 - 0.2e1 * t2 * t19; 0, 0, 0, 0, 0, t73, t72, t45, t29, -t30, t44 * t78 - t58, -t42 * t78 - t11, (-t27 * t42 - t28 * t44) * pkin(2), (t11 * t42 - t44 * t58) * pkin(2), t20 * t47, t20 * t50 - t67, t66, t26, 0, t37 * t19 - t27 * t65 - t8 * t50, t37 * t20 - t27 * t60 + t8 * t47, t13 * t63, (-t12 * t49 - t77) * t47, t53 - t76, -t54 + t61, -t19 * t50, -t1 * t50 + t17 * t19 + (t12 * t36 + t82) * t47, -t18 * t19 + t2 * t50 + (t13 * t36 + t81) * t47; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t79, -0.2e1 * t80, 0, (t42 ^ 2 + t44 ^ 2) * pkin(2) ^ 2, t40, t56, 0, 0, 0, -0.2e1 * t37 * t50, t37 * t85, t41 * t40, -0.2e1 * t40 * t69, -0.2e1 * t47 * t62, t46 * t56, t50 ^ 2, -0.2e1 * t17 * t50 + 0.2e1 * t40 * t75, 0.2e1 * t40 * t36 * t49 + 0.2e1 * t18 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t28, 0, t32, 0, 0, 0, 0, 0, t26, -t66, 0, 0, 0, 0, 0, -t54 - t61, -t53 - t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t27, t6, -t7, t77, -t46 * t12 + t13 * t49, t71, t64, 0, -pkin(4) * t12 - pkin(9) * t71 - t81, -pkin(4) * t13 - pkin(9) * t64 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t50, 0, -t65, -t60, t46 * t63, (-t39 + t41) * t47, -t68, -t62, 0, -t36 * t63 + (-pkin(4) * t47 + pkin(9) * t50) * t46, pkin(9) * t62 + (t75 - t83) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t47, 0, 0, 0, 0, 0, t62, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t39, 0.2e1 * t69, 0, 0, 0, 0.2e1 * t83, -0.2e1 * pkin(4) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, t19, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t70, -t50, t17, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t49, 0, -t46 * pkin(9), -t49 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t10;
