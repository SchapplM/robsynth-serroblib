% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:04
% EndTime: 2019-12-31 20:14:06
% DurationCPUTime: 0.68s
% Computational Cost: add. (316->79), mult. (579->128), div. (0->0), fcn. (520->4), ass. (0->60)
t43 = sin(qJ(4));
t38 = t43 ^ 2;
t45 = cos(qJ(4));
t40 = t45 ^ 2;
t24 = t38 + t40;
t44 = sin(qJ(2));
t39 = t44 ^ 2;
t46 = cos(qJ(2));
t41 = t46 ^ 2;
t70 = t39 + t41;
t69 = -0.2e1 * t44;
t68 = 0.2e1 * t46;
t67 = 2 * qJ(3);
t47 = -pkin(2) - pkin(7);
t66 = t44 * pkin(4);
t53 = -t44 * qJ(3) - pkin(1);
t10 = t46 * t47 + t53;
t34 = t44 * pkin(6);
t20 = pkin(3) * t44 + t34;
t6 = t45 * t10 + t43 * t20;
t65 = t43 * t46;
t64 = t43 * t47;
t63 = t44 * t46;
t62 = t44 * t47;
t61 = t45 * t43;
t60 = t45 * t46;
t32 = t45 * t47;
t59 = t24 * t47 ^ 2;
t58 = t70 * pkin(6) ^ 2;
t36 = t46 * pkin(6);
t21 = t46 * pkin(3) + t36;
t57 = t44 * qJ(5);
t56 = t46 * qJ(3);
t55 = t41 * t61;
t54 = t44 * t60;
t52 = t10 * t43 - t45 * t20;
t3 = t57 + t6;
t4 = t52 - t66;
t1 = t3 * t43 - t4 * t45;
t2 = t43 * t6 - t45 * t52;
t51 = -pkin(2) * t44 + t56;
t18 = pkin(4) * t45 + qJ(5) * t43;
t50 = pkin(4) * t43 - qJ(5) * t45;
t48 = qJ(3) ^ 2;
t29 = t45 * t44;
t28 = t40 * t41;
t27 = t43 * t44;
t26 = t38 * t41;
t25 = 0.2e1 * t63;
t23 = t44 * t32;
t22 = t43 * t60;
t19 = -0.2e1 * t43 * t63;
t17 = -pkin(2) * t46 + t53;
t16 = qJ(3) + t50;
t15 = 0.2e1 * t70 * pkin(6);
t14 = t24 * t47;
t11 = (-t38 + t40) * t46;
t9 = 0.2e1 * t14;
t7 = t18 * t46 + t21;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t39, t25, 0, t41, 0, 0, pkin(1) * t68, pkin(1) * t69, t15, pkin(1) ^ 2 + t58, 0, 0, 0, t39, t25, t41, t15, t17 * t68, t17 * t69, t17 ^ 2 + t58, t26, 0.2e1 * t55, t19, t28, -0.2e1 * t54, t39, 0.2e1 * t21 * t60 - 0.2e1 * t44 * t52, -0.2e1 * t21 * t65 - 0.2e1 * t44 * t6, (-t43 * t52 - t45 * t6) * t68, t21 ^ 2 + t52 ^ 2 + t6 ^ 2, t26, t19, -0.2e1 * t55, t39, 0.2e1 * t54, t28, -0.2e1 * t4 * t44 + 0.2e1 * t60 * t7, (-t3 * t45 - t4 * t43) * t68, 0.2e1 * t3 * t44 + 0.2e1 * t65 * t7, t3 ^ 2 + t4 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t46, 0, -t34, -t36, 0, 0, 0, -t44, -t46, 0, 0, 0, t51, t34, t36, t51 * pkin(6), -t22, -t11, t29, t22, -t27, 0, t21 * t43 + t45 * t56 + t23, t21 * t45 + (-t56 - t62) * t43, -t2, t21 * qJ(3) + t2 * t47, -t22, t29, t11, 0, t27, t22, t16 * t60 + t43 * t7 + t23, -t1, -t7 * t45 + (t16 * t46 + t62) * t43, t1 * t47 + t7 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(2), t67, pkin(2) ^ 2 + t48, t40, -0.2e1 * t61, 0, t38, 0, 0, t43 * t67, t45 * t67, -t9, t48 + t59, t40, 0, 0.2e1 * t61, 0, 0, t38, 0.2e1 * t16 * t43, -t9, -0.2e1 * t16 * t45, t16 ^ 2 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, t34, 0, 0, 0, 0, 0, 0, t29, -t27, 0, t2, 0, 0, 0, 0, 0, 0, t29, 0, t27, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t24, t14, 0, 0, 0, 0, 0, 0, 0, -t24, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, 0, -t60, t44, -t52, -t6, 0, 0, 0, -t65, 0, t44, t60, 0, -t52 + 0.2e1 * t66, t50 * t46, 0.2e1 * t57 + t6, -pkin(4) * t4 + qJ(5) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, -t43, 0, t32, -t64, 0, 0, 0, t45, 0, 0, t43, 0, t32, -t18, t64, t18 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t43, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t43, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t65, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t5;
