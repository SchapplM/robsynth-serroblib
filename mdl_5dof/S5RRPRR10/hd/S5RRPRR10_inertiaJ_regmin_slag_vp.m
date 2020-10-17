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
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 20:26:39
% EndTime: 2019-12-31 20:26:41
% DurationCPUTime: 0.48s
% Computational Cost: add. (653->93), mult. (1716->203), div. (0->0), fcn. (1964->10), ass. (0->75)
t41 = sin(pkin(10));
t42 = sin(pkin(5));
t43 = cos(pkin(10));
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t27 = (t41 * t50 + t43 * t47) * t42;
t44 = cos(pkin(5));
t46 = sin(qJ(4));
t49 = cos(qJ(4));
t18 = t27 * t46 - t44 * t49;
t81 = -0.2e1 * t18;
t80 = 0.2e1 * t46;
t79 = pkin(1) * t47;
t48 = cos(qJ(5));
t78 = pkin(4) * t48;
t70 = t42 * t50;
t71 = t42 * t47;
t26 = t41 * t71 - t43 * t70;
t31 = (-pkin(2) * t50 - pkin(1)) * t42;
t14 = t26 * pkin(3) - t27 * pkin(8) + t31;
t33 = t44 * t50 * pkin(1);
t57 = pkin(7) + qJ(3);
t20 = t44 * pkin(2) - t57 * t71 + t33;
t54 = t44 * t79;
t23 = t57 * t70 + t54;
t11 = t41 * t20 + t43 * t23;
t9 = t44 * pkin(8) + t11;
t6 = t49 * t14 - t46 * t9;
t3 = -t26 * pkin(4) - t6;
t45 = sin(qJ(5));
t77 = t3 * t45;
t76 = t3 * t48;
t19 = t27 * t49 + t44 * t46;
t13 = t19 * t48 + t26 * t45;
t75 = t13 * t45;
t74 = t13 * t49;
t35 = t41 * pkin(2) + pkin(8);
t73 = t35 * t45;
t37 = t42 ^ 2;
t72 = t37 * t50;
t69 = t45 * t18;
t68 = t45 * t46;
t67 = t45 * t48;
t66 = t45 * t49;
t65 = t46 * t18;
t64 = t46 * t26;
t63 = t46 * t35;
t62 = t48 * t18;
t61 = t48 * t46;
t60 = t48 * t49;
t12 = t19 * t45 - t26 * t48;
t59 = t49 * t12;
t58 = t49 * t35;
t56 = 0.2e1 * t42 * t44;
t55 = t49 * t80;
t53 = t45 * t65;
t52 = t18 * t61;
t36 = -t43 * pkin(2) - pkin(3);
t10 = t43 * t20 - t41 * t23;
t7 = t46 * t14 + t49 * t9;
t8 = -t44 * pkin(3) - t10;
t40 = t48 ^ 2;
t39 = t46 ^ 2;
t38 = t45 ^ 2;
t30 = -t49 * pkin(4) - t46 * pkin(9) + t36;
t29 = pkin(7) * t70 + t54;
t28 = -pkin(7) * t71 + t33;
t24 = t49 * t26;
t17 = t45 * t30 + t48 * t58;
t16 = t48 * t30 - t45 * t58;
t5 = t18 * pkin(4) - t19 * pkin(9) + t8;
t4 = t26 * pkin(9) + t7;
t2 = t48 * t4 + t45 * t5;
t1 = -t45 * t4 + t48 * t5;
t15 = [1, 0, 0, t37 * t47 ^ 2, 0.2e1 * t47 * t72, t47 * t56, t50 * t56, t44 ^ 2, 0.2e1 * pkin(1) * t72 + 0.2e1 * t28 * t44, -0.2e1 * t29 * t44 - 0.2e1 * t37 * t79, -0.2e1 * t10 * t27 - 0.2e1 * t11 * t26, t10 ^ 2 + t11 ^ 2 + t31 ^ 2, t19 ^ 2, t19 * t81, 0.2e1 * t19 * t26, t26 * t81, t26 ^ 2, 0.2e1 * t8 * t18 + 0.2e1 * t6 * t26, 0.2e1 * t8 * t19 - 0.2e1 * t7 * t26, t13 ^ 2, -0.2e1 * t13 * t12, 0.2e1 * t13 * t18, t12 * t81, t18 ^ 2, 0.2e1 * t1 * t18 + 0.2e1 * t3 * t12, 0.2e1 * t3 * t13 - 0.2e1 * t2 * t18; 0, 0, 0, 0, 0, t71, t70, t44, t28, -t29, (-t26 * t41 - t27 * t43) * pkin(2), (t10 * t43 + t11 * t41) * pkin(2), t19 * t46, t19 * t49 - t65, t64, t24, 0, t36 * t18 - t26 * t63 - t8 * t49, t36 * t19 - t26 * t58 + t8 * t46, t13 * t61, (-t12 * t48 - t75) * t46, t52 - t74, -t53 + t59, -t18 * t49, -t1 * t49 + t16 * t18 + (t12 * t35 + t77) * t46, -t17 * t18 + t2 * t49 + (t13 * t35 + t76) * t46; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t41 ^ 2 + t43 ^ 2) * pkin(2) ^ 2, t39, t55, 0, 0, 0, -0.2e1 * t36 * t49, t36 * t80, t40 * t39, -0.2e1 * t39 * t67, -0.2e1 * t46 * t60, t45 * t55, t49 ^ 2, -0.2e1 * t16 * t49 + 0.2e1 * t39 * t73, 0.2e1 * t39 * t35 * t48 + 0.2e1 * t17 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, t24, -t64, 0, 0, 0, 0, 0, -t53 - t59, -t52 - t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t26, t6, -t7, t75, -t45 * t12 + t13 * t48, t69, t62, 0, -pkin(4) * t12 - pkin(9) * t69 - t76, -pkin(4) * t13 - pkin(9) * t62 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t49, 0, -t63, -t58, t45 * t61, (-t38 + t40) * t46, -t66, -t60, 0, -t35 * t61 + (-pkin(4) * t46 + pkin(9) * t49) * t45, pkin(9) * t60 + (t73 - t78) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t46, 0, 0, 0, 0, 0, t60, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t38, 0.2e1 * t67, 0, 0, 0, 0.2e1 * t78, -0.2e1 * pkin(4) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, t18, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t68, -t49, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t48, 0, -t45 * pkin(9), -t48 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t15;
