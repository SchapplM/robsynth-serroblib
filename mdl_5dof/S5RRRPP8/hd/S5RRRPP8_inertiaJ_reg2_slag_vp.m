% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:24
% EndTime: 2019-12-31 21:09:27
% DurationCPUTime: 0.83s
% Computational Cost: add. (354->105), mult. (745->178), div. (0->0), fcn. (675->4), ass. (0->69)
t46 = sin(qJ(3));
t40 = t46 ^ 2;
t48 = cos(qJ(3));
t42 = t48 ^ 2;
t80 = t40 + t42;
t32 = qJ(4) * t48;
t79 = -pkin(3) * t46 + t32;
t44 = pkin(3) + qJ(5);
t64 = t46 * qJ(4);
t78 = -t44 * t48 - t64;
t77 = -0.2e1 * t46;
t47 = sin(qJ(2));
t76 = 0.2e1 * t47;
t75 = pkin(2) * t48;
t49 = cos(qJ(2));
t73 = pkin(7) * t49;
t41 = t47 ^ 2;
t72 = t41 * pkin(6);
t36 = t47 * pkin(6);
t27 = t46 * t47;
t70 = t46 * t48;
t69 = t46 * t49;
t68 = t47 * t49;
t29 = t48 * t47;
t30 = t48 * t49;
t67 = pkin(3) * t27 + t36;
t15 = -t49 * pkin(2) - t47 * pkin(7) - pkin(1);
t66 = pkin(6) * t69 - t48 * t15;
t7 = pkin(6) * t30 + t46 * t15;
t65 = t80 * pkin(7) ^ 2;
t63 = t49 * qJ(4);
t62 = t46 * t68;
t61 = t41 * t70;
t60 = t47 * t30;
t39 = t49 * pkin(3);
t5 = t39 + t66;
t4 = t63 - t7;
t59 = -t4 * t48 + t5 * t46;
t58 = t46 * t66 + t7 * t48;
t56 = -t48 * pkin(3) - t64;
t14 = -pkin(2) + t56;
t57 = -t14 * t47 - t73;
t55 = -pkin(4) * t27 + t7;
t54 = -pkin(4) * t29 - t5;
t53 = pkin(6) ^ 2;
t51 = qJ(4) ^ 2;
t50 = 0.2e1 * qJ(4);
t43 = t49 ^ 2;
t38 = t48 * pkin(7);
t35 = t41 * t53;
t34 = t46 * pkin(7);
t31 = -0.2e1 * t63;
t28 = t42 * t41;
t26 = t40 * t41;
t25 = 0.2e1 * t70;
t21 = t46 * t29;
t20 = t48 * pkin(4) + t38;
t19 = t46 * pkin(4) + t34;
t18 = -0.2e1 * t60;
t17 = -0.2e1 * t61;
t16 = -0.2e1 * t62;
t13 = 0.2e1 * t80 * pkin(7);
t12 = (-t40 + t42) * t47;
t9 = -pkin(2) + t78;
t8 = -t47 * t32 + t67;
t3 = (qJ(5) * t46 - t32) * t47 + t67;
t2 = t55 - t63;
t1 = t49 * qJ(5) - t54;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t41, 0.2e1 * t68, 0, t43, 0, 0, 0.2e1 * pkin(1) * t49, -0.2e1 * pkin(1) * t47, 0.2e1 * (t41 + t43) * pkin(6), pkin(1) ^ 2 + t43 * t53 + t35, t28, t17, t18, t26, 0.2e1 * t62, t43, 0.2e1 * t46 * t72 + 0.2e1 * t49 * t66, 0.2e1 * t48 * t72 + 0.2e1 * t7 * t49, (-t46 * t7 + t48 * t66) * t76, t66 ^ 2 + t7 ^ 2 + t35, t43, 0.2e1 * t60, t16, t28, t17, t26, (t4 * t46 + t48 * t5) * t76, -0.2e1 * t8 * t27 - 0.2e1 * t5 * t49, -0.2e1 * t8 * t29 + 0.2e1 * t4 * t49, t4 ^ 2 + t5 ^ 2 + t8 ^ 2, t43, t16, t18, t26, 0.2e1 * t61, t28, (t1 * t48 - t2 * t46) * t76, -0.2e1 * t2 * t49 - 0.2e1 * t3 * t29, 0.2e1 * t1 * t49 + 0.2e1 * t3 * t27, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t49, 0, -t36, -t49 * pkin(6), 0, 0, t21, t12, -t69, -t21, -t30, 0, -pkin(6) * t29 + (-pkin(2) * t47 + t73) * t46, pkin(7) * t30 + (pkin(6) * t46 - t75) * t47, t58, -pkin(2) * t36 + t58 * pkin(7), 0, t69, t30, t21, t12, -t21, t59, t57 * t46 + t8 * t48, -t8 * t46 + t57 * t48, t59 * pkin(7) + t8 * t14, 0, t30, -t69, -t21, -t12, t21, (t19 * t47 + t2) * t48 + (-t20 * t47 + t1) * t46, -t20 * t49 - t9 * t29 - t3 * t46, t19 * t49 + t9 * t27 - t3 * t48, t1 * t19 + t2 * t20 + t3 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t40, t25, 0, t42, 0, 0, 0.2e1 * t75, pkin(2) * t77, t13, pkin(2) ^ 2 + t65, 0, 0, 0, t40, t25, t42, t13, 0.2e1 * t14 * t48, t14 * t77, t14 ^ 2 + t65, 0, 0, 0, t42, -0.2e1 * t70, t40, 0.2e1 * t19 * t46 + 0.2e1 * t20 * t48, t9 * t77, -0.2e1 * t9 * t48, t19 ^ 2 + t20 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t27, -t49, -t66, -t7, 0, 0, -t49, -t29, t27, 0, 0, 0, t56 * t47, 0.2e1 * t39 + t66, t31 + t7, -t5 * pkin(3) - t4 * qJ(4), -t49, t27, t29, 0, 0, 0, t78 * t47, t31 + t55, (-qJ(5) - t44) * t49 + t54, t2 * qJ(4) - t1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, t48, 0, -t34, -t38, 0, 0, 0, -t46, -t48, 0, 0, 0, t79, t34, t38, t79 * pkin(7), 0, -t48, t46, 0, 0, 0, -t44 * t46 + t32, t20, -t19, t20 * qJ(4) - t19 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(3), t50, pkin(3) ^ 2 + t51, 1, 0, 0, 0, 0, 0, 0, t50, 0.2e1 * t44, t44 ^ 2 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t49, 0, t5, 0, 0, 0, 0, 0, 0, t29, 0, t49, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, t34, 0, 0, 0, 0, 0, 0, t46, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -1, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t49, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
