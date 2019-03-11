% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t42 = sin(qJ(4));
t48 = pkin(4) + pkin(5);
t45 = cos(qJ(4));
t63 = t45 * qJ(5);
t87 = t48 * t42 - t63;
t64 = t42 * qJ(5);
t86 = t45 * t48 + t64;
t22 = pkin(3) + t86;
t85 = 0.2e1 * t22;
t56 = -t45 * pkin(4) - t64;
t25 = -pkin(3) + t56;
t84 = -0.2e1 * t25;
t43 = sin(qJ(3));
t83 = -0.2e1 * t43;
t82 = 0.2e1 * t43;
t46 = cos(qJ(3));
t81 = 0.2e1 * t46;
t80 = pkin(3) * t42;
t79 = pkin(3) * t45;
t78 = pkin(8) * t42;
t77 = pkin(8) * t45;
t76 = t42 * pkin(9);
t60 = cos(pkin(6));
t41 = sin(pkin(6));
t72 = t41 * sin(qJ(2));
t21 = t60 * t43 + t46 * t72;
t71 = t41 * cos(qJ(2));
t10 = t21 * t42 + t45 * t71;
t19 = t43 * t72 - t60 * t46;
t70 = t42 * t43;
t75 = t10 * t46 + t19 * t70;
t74 = t19 * t42;
t73 = t19 * t45;
t69 = t42 * t45;
t68 = t42 * t46;
t34 = t45 * t43;
t67 = t45 * t46;
t26 = -t46 * pkin(3) - t43 * pkin(9) - pkin(2);
t66 = pkin(8) * t68 - t45 * t26;
t17 = pkin(8) * t67 + t42 * t26;
t38 = t42 ^ 2;
t40 = t45 ^ 2;
t65 = t38 + t40;
t62 = t45 * qJ(6);
t61 = t46 * qJ(5);
t59 = t43 * t81;
t37 = t46 * pkin(4);
t14 = t37 + t66;
t58 = -0.2e1 * t61 + t17;
t11 = t21 * t45 - t42 * t71;
t57 = t10 ^ 2 + t11 ^ 2 + t19 ^ 2;
t13 = -t61 + t17;
t55 = -pkin(4) * t42 + t63;
t2 = t10 * t42 + t11 * t45;
t54 = t13 * t45 + t14 * t42;
t3 = t11 * t46 + t19 * t34;
t53 = t43 * t62 - t14;
t51 = qJ(5) ^ 2;
t50 = 0.2e1 * qJ(5);
t39 = t43 ^ 2;
t36 = t45 * pkin(9);
t31 = pkin(9) * t68;
t29 = qJ(6) * t70;
t28 = t36 - t62;
t27 = (pkin(9) - qJ(6)) * t42;
t18 = (pkin(8) - t55) * t43;
t12 = (-pkin(8) - t87) * t43;
t7 = t11 * qJ(5);
t5 = t13 + t29;
t4 = t46 * pkin(5) - t53;
t1 = (-t10 * t45 + t11 * t42) * t43;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, t57; 0, 0, t71, -t72, 0, 0, 0, 0, 0, t46 * t71, -t43 * t71, 0, 0, 0, 0, 0, t75, t3, t75, -t1, -t3, t10 * t14 + t11 * t13 + t19 * t18, t75, -t3, t1, t10 * t4 + t11 * t5 - t19 * t12; 0, 1, 0, 0, t39, t59, 0, 0, 0, pkin(2) * t81, pkin(2) * t83, t40 * t39, -0.2e1 * t39 * t69, t67 * t83, t42 * t59, t46 ^ 2, 0.2e1 * t39 * t78 + 0.2e1 * t46 * t66, 0.2e1 * t17 * t46 + 0.2e1 * t39 * t77, 0.2e1 * t14 * t46 + 0.2e1 * t18 * t70 (-t13 * t42 + t14 * t45) * t82, -0.2e1 * t13 * t46 - 0.2e1 * t18 * t34, t13 ^ 2 + t14 ^ 2 + t18 ^ 2, -0.2e1 * t12 * t70 + 0.2e1 * t4 * t46, 0.2e1 * t12 * t34 - 0.2e1 * t5 * t46 (-t4 * t45 + t42 * t5) * t82, t12 ^ 2 + t4 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t21, 0, 0, 0, 0, 0, -t73, t74, -t73, t2, -t74, pkin(9) * t2 + t19 * t25, -t73, -t74, -t2, t10 * t27 + t11 * t28 - t19 * t22; 0, 0, 0, 0, 0, 0, t43, t46, 0, -t43 * pkin(8), -t46 * pkin(8), t42 * t34 (-t38 + t40) * t43, -t68, -t67, 0, t31 + (-t77 - t80) * t43, pkin(9) * t67 + (t78 - t79) * t43, -t18 * t45 + t25 * t70 + t31, t54, -t18 * t42 + (-pkin(9) * t46 - t25 * t43) * t45, pkin(9) * t54 + t18 * t25, t12 * t45 - t22 * t70 + t27 * t46, t12 * t42 + t22 * t34 - t28 * t46 (-t27 * t43 - t5) * t45 + (t28 * t43 - t4) * t42, t12 * t22 + t4 * t27 + t5 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t38, 0.2e1 * t69, 0, 0, 0, 0.2e1 * t79, -0.2e1 * t80, t45 * t84, 0.2e1 * t65 * pkin(9), t42 * t84, t65 * pkin(9) ^ 2 + t25 ^ 2, t45 * t85, t42 * t85, -0.2e1 * t27 * t42 - 0.2e1 * t28 * t45, t22 ^ 2 + t27 ^ 2 + t28 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, -t10, 0, t11, -t10 * pkin(4) + t7, -t10, t11, 0, -t10 * t48 + t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t70, -t46, -t66, -t17, -0.2e1 * t37 - t66, t56 * t43, t58, -t14 * pkin(4) + t13 * qJ(5) (-pkin(5) - t48) * t46 + t53, t29 + t58, t86 * t43, t5 * qJ(5) - t4 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t45, 0, -t76, -t36, -t76, t55, t36, t55 * pkin(9), -t27, t28, t87, t28 * qJ(5) - t27 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, t50, pkin(4) ^ 2 + t51, 0.2e1 * t48, t50, 0, t48 ^ 2 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t34, 0, t14, t46, 0, -t34, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t76, 0, 0, -t42, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), -1, 0, 0, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t34, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t42, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
