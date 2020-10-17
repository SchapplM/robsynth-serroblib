% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR13_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:33
% EndTime: 2019-12-31 21:46:35
% DurationCPUTime: 0.62s
% Computational Cost: add. (523->107), mult. (1293->220), div. (0->0), fcn. (1405->8), ass. (0->83)
t44 = cos(pkin(5));
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t43 = sin(pkin(5));
t47 = sin(qJ(2));
t77 = t43 * t47;
t25 = t44 * t46 + t49 * t77;
t88 = -0.2e1 * t25;
t87 = 0.2e1 * t25;
t86 = -0.2e1 * t46;
t85 = 0.2e1 * t49;
t84 = 2 * qJ(4);
t83 = pkin(3) + pkin(9);
t82 = pkin(1) * t47;
t50 = cos(qJ(2));
t81 = pkin(1) * t50;
t24 = -t44 * t49 + t46 * t77;
t45 = sin(qJ(5));
t48 = cos(qJ(5));
t76 = t43 * t50;
t14 = t24 * t45 - t48 * t76;
t80 = t14 * t48;
t79 = t25 * t45;
t21 = t25 * t46;
t38 = t43 ^ 2;
t78 = t38 * t50;
t75 = t44 * t47;
t74 = t45 * t46;
t73 = t45 * t49;
t72 = t45 * t83;
t71 = t46 * t49;
t70 = t48 * t45;
t69 = t48 * t49;
t68 = t48 * t83;
t63 = pkin(7) * t76;
t19 = t63 + (pkin(8) + t82) * t44;
t20 = (-pkin(2) * t50 - pkin(8) * t47 - pkin(1)) * t43;
t10 = t49 * t19 + t46 * t20;
t40 = t46 ^ 2;
t42 = t49 ^ 2;
t67 = t40 + t42;
t66 = qJ(4) * t49;
t65 = 0.2e1 * t76;
t64 = -0.2e1 * t71;
t62 = t46 * t76;
t61 = t49 * t76;
t60 = qJ(4) * t76;
t59 = -t46 * qJ(4) - pkin(2);
t9 = -t46 * t19 + t49 * t20;
t58 = pkin(8) * t62;
t57 = pkin(8) * t61;
t33 = pkin(3) * t76;
t8 = t33 - t9;
t7 = t60 - t10;
t56 = t8 * t46 - t7 * t49;
t55 = -pkin(3) * t46 + t66;
t54 = -t46 * t83 + t66;
t32 = pkin(7) * t77;
t18 = t32 + (-pkin(2) - t81) * t44;
t53 = -t25 * qJ(4) + t18;
t41 = t48 ^ 2;
t39 = t45 ^ 2;
t37 = t49 * pkin(8);
t36 = t46 * pkin(8);
t35 = t48 * t46;
t31 = t49 * pkin(4) + t37;
t30 = t46 * pkin(4) + t36;
t29 = -t49 * pkin(3) + t59;
t28 = -t83 * t49 + t59;
t27 = pkin(1) * t75 + t63;
t26 = t44 * t81 - t32;
t23 = t25 ^ 2;
t22 = t25 * t48;
t13 = t24 * t48 + t45 * t76;
t12 = t48 * t28 + t45 * t30;
t11 = -t45 * t28 + t48 * t30;
t6 = t24 * pkin(3) + t53;
t5 = -t24 * pkin(4) - t7;
t4 = t83 * t24 + t53;
t3 = t25 * pkin(4) + pkin(9) * t76 + t8;
t2 = t45 * t3 + t48 * t4;
t1 = t48 * t3 - t45 * t4;
t15 = [1, 0, 0, t38 * t47 ^ 2, 0.2e1 * t47 * t78, 0.2e1 * t43 * t75, t44 * t65, t44 ^ 2, 0.2e1 * pkin(1) * t78 + 0.2e1 * t26 * t44, -0.2e1 * t27 * t44 - 0.2e1 * t38 * t82, t23, t24 * t88, t76 * t88, t24 * t65, t38 * t50 ^ 2, 0.2e1 * t18 * t24 - 0.2e1 * t9 * t76, 0.2e1 * t10 * t76 + 0.2e1 * t18 * t25, 0.2e1 * t7 * t24 + 0.2e1 * t8 * t25, -0.2e1 * t6 * t24 - 0.2e1 * t8 * t76, -0.2e1 * t6 * t25 + 0.2e1 * t7 * t76, t6 ^ 2 + t7 ^ 2 + t8 ^ 2, t14 ^ 2, 0.2e1 * t14 * t13, t14 * t87, t13 * t87, t23, 0.2e1 * t1 * t25 - 0.2e1 * t5 * t13, 0.2e1 * t5 * t14 - 0.2e1 * t2 * t25; 0, 0, 0, 0, 0, t77, t76, t44, t26, -t27, t21, -t46 * t24 + t25 * t49, -t62, -t61, 0, -pkin(2) * t24 - t18 * t49 + t58, -pkin(2) * t25 + t18 * t46 + t57, (-t24 * t49 + t21) * pkin(8) + t56, -t29 * t24 + t6 * t49 - t58, -t29 * t25 - t6 * t46 - t57, t56 * pkin(8) + t6 * t29, -t14 * t73, (-t13 * t45 - t80) * t49, t14 * t46 - t25 * t73, t13 * t46 - t25 * t69, t21, t1 * t46 + t11 * t25 - t31 * t13 + t5 * t69, -t12 * t25 + t31 * t14 - t2 * t46 - t5 * t73; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t40, 0.2e1 * t71, 0, 0, 0, pkin(2) * t85, pkin(2) * t86, 0.2e1 * t67 * pkin(8), t29 * t85, t29 * t86, t67 * pkin(8) ^ 2 + t29 ^ 2, t39 * t42, 0.2e1 * t42 * t70, t45 * t64, t48 * t64, t40, 0.2e1 * t11 * t46 + 0.2e1 * t31 * t69, -0.2e1 * t12 * t46 - 0.2e1 * t31 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, -t76, t9, -t10, -t25 * pkin(3) - qJ(4) * t24, 0.2e1 * t33 - t9, -0.2e1 * t60 + t10, -t8 * pkin(3) - t7 * qJ(4), t80, t48 * t13 - t14 * t45, t22, -t79, 0, -qJ(4) * t13 - t25 * t68 + t5 * t45, qJ(4) * t14 + t25 * t72 + t5 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t49, 0, -t36, -t37, t55, t36, t37, t55 * pkin(8), -t45 * t69, (t39 - t41) * t49, t35, -t74, 0, t31 * t45 + t54 * t48, t31 * t48 - t54 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t84, pkin(3) ^ 2 + (qJ(4) ^ 2), t41, -0.2e1 * t70, 0, 0, 0, t45 * t84, t48 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t76, 0, t8, 0, 0, 0, 0, 0, t22, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, t36, 0, 0, 0, 0, 0, t35, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t13, t25, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t69, t46, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t45, 0, -t68, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t15;
