% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_inertiaJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:51:00
% EndTime: 2019-12-05 18:51:03
% DurationCPUTime: 0.82s
% Computational Cost: add. (492->67), mult. (1038->150), div. (0->0), fcn. (1223->8), ass. (0->58)
t48 = cos(qJ(3));
t37 = t48 * pkin(2);
t32 = t37 + pkin(3);
t43 = sin(qJ(4));
t47 = cos(qJ(4));
t44 = sin(qJ(3));
t69 = t44 * pkin(2);
t58 = t47 * t69;
t19 = t43 * t32 + t58;
t17 = pkin(6) + t19;
t42 = sin(qJ(5));
t40 = t42 ^ 2;
t46 = cos(qJ(5));
t41 = t46 ^ 2;
t61 = t40 + t41;
t64 = t61 * t17;
t70 = t43 * pkin(3);
t30 = pkin(6) + t70;
t77 = t61 * t30;
t45 = sin(qJ(2));
t49 = cos(qJ(2));
t21 = t44 * t45 - t48 * t49;
t22 = t44 * t49 + t48 * t45;
t10 = t43 * t21 - t47 * t22;
t76 = -0.2e1 * t10;
t8 = -t47 * t21 - t43 * t22;
t75 = t8 ^ 2;
t33 = t49 * pkin(2) + pkin(1);
t14 = -t21 * pkin(3) + t33;
t74 = 0.2e1 * t14;
t73 = -0.2e1 * t33;
t72 = 0.2e1 * t49;
t71 = pkin(4) * t42;
t5 = t42 * t8;
t6 = t46 * t8;
t56 = -t47 * t32 + t43 * t69;
t16 = -pkin(4) + t56;
t68 = t16 * t46;
t36 = t47 * pkin(3);
t31 = -t36 - pkin(4);
t67 = t31 * t46;
t66 = t42 * t46;
t65 = t46 * t10;
t62 = pkin(6) * t61;
t60 = -0.2e1 * t5;
t59 = 0.2e1 * t6;
t55 = -pkin(4) * t10 - pkin(6) * t8;
t54 = t10 * t16 - t17 * t8;
t53 = t10 * t31 - t30 * t8;
t39 = pkin(4) * t46;
t28 = 0.2e1 * t66;
t27 = t31 * t42;
t13 = t16 * t42;
t7 = t10 ^ 2;
t4 = t42 * t65;
t3 = (-t40 + t41) * t10;
t2 = t8 * pkin(4) - t10 * pkin(6) + t14;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t45 ^ 2, t45 * t72, 0, t49 ^ 2, 0, 0, pkin(1) * t72, -0.2e1 * pkin(1) * t45, 0, pkin(1) ^ 2, t22 ^ 2, -0.2e1 * t21 * t22, 0, t21 ^ 2, 0, 0, t21 * t73, t22 * t73, 0, t33 ^ 2, t7, t8 * t76, 0, t75, 0, 0, t8 * t74, t10 * t74, 0, t14 ^ 2, t41 * t7, -0.2e1 * t7 * t66, t10 * t59, t40 * t7, t10 * t60, t75, t2 * t59, t2 * t60, t61 * t2 * t76, t61 * t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, 0, -t49, 0, 0, 0, 0, 0, 0, 0, -t22, 0, t21, 0, 0, 0, (t21 * t44 + t22 * t48) * pkin(2), 0, 0, 0, t10, 0, -t8, 0, 0, 0, t10 * t56 - t19 * t8, 0, t4, t3, t5, -t4, t6, 0, t54 * t42, t54 * t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t69, 0, (t44 ^ 2 + t48 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t56, -0.2e1 * t19, 0, t19 ^ 2 + t56 ^ 2, t40, t28, 0, t41, 0, 0, -0.2e1 * t68, 0.2e1 * t13, 0.2e1 * t64, t61 * t17 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, t21, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t8, 0, 0, 0, (-t10 * t47 - t43 * t8) * pkin(3), 0, t4, t3, t5, -t4, t6, 0, t53 * t42, t53 * t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t37, -t69, 0, 0, 0, 0, 0, 0, 0, 1, t36 - t56, -t58 + (-pkin(3) - t32) * t43, 0, (t19 * t43 - t47 * t56) * pkin(3), t40, t28, 0, t41, 0, 0, (-t16 - t31) * t46, t27 + t13, t77 + t64, t16 * t31 + t17 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t36, -0.2e1 * t70, 0, (t43 ^ 2 + t47 ^ 2) * pkin(3) ^ 2, t40, t28, 0, t41, 0, 0, -0.2e1 * t67, 0.2e1 * t27, 0.2e1 * t77, t61 * t30 ^ 2 + t31 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t8, 0, 0, 0, 0, 0, t4, t3, t5, -t4, t6, 0, t55 * t42, t55 * t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t56, -t19, 0, 0, t40, t28, 0, t41, 0, 0, t39 - t68, t13 - t71, t62 + t64, -t16 * pkin(4) + pkin(6) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t36, -t70, 0, 0, t40, t28, 0, t41, 0, 0, t39 - t67, t27 - t71, t62 + t77, -t31 * pkin(4) + pkin(6) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t40, t28, 0, t41, 0, 0, 0.2e1 * t39, -0.2e1 * t71, 0.2e1 * t62, t61 * pkin(6) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, -t42 * t10, t8, t46 * t2, -t42 * t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t46, 0, -t42 * t17, -t46 * t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t46, 0, -t42 * t30, -t46 * t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t46, 0, -t42 * pkin(6), -t46 * pkin(6), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
