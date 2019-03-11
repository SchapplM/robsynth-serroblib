% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t66 = sin(qJ(2));
t50 = (-qJ(3) - pkin(7)) * t66;
t80 = cos(qJ(2));
t73 = t80 * pkin(7);
t51 = qJ(3) * t80 + t73;
t61 = sin(pkin(10));
t63 = cos(pkin(10));
t34 = -t63 * t50 + t61 * t51;
t90 = t34 ^ 2;
t60 = sin(pkin(11));
t62 = cos(pkin(11));
t65 = sin(qJ(5));
t68 = cos(qJ(5));
t45 = t65 * t60 - t68 * t62;
t56 = -t63 * pkin(2) - pkin(3);
t49 = -t62 * pkin(4) + t56;
t37 = t45 * pkin(5) + t49;
t89 = 0.2e1 * t37;
t44 = t61 * t66 - t63 * t80;
t88 = -0.2e1 * t44;
t87 = 0.2e1 * t44;
t86 = 0.2e1 * t49;
t85 = t44 * pkin(5);
t64 = sin(qJ(6));
t84 = t64 * pkin(5);
t67 = cos(qJ(6));
t83 = t67 * pkin(5);
t46 = t61 * t80 + t63 * t66;
t47 = t68 * t60 + t65 * t62;
t23 = t47 * t46;
t57 = -pkin(2) * t80 - pkin(1);
t29 = t44 * pkin(3) - t46 * qJ(4) + t57;
t36 = t61 * t50 + t63 * t51;
t16 = t62 * t29 - t60 * t36;
t76 = t62 * t46;
t11 = t44 * pkin(4) - pkin(8) * t76 + t16;
t17 = t60 * t29 + t62 * t36;
t77 = t60 * t46;
t12 = -pkin(8) * t77 + t17;
t7 = t65 * t11 + t68 * t12;
t5 = -t23 * pkin(9) + t7;
t82 = t67 * t5;
t53 = t61 * pkin(2) + qJ(4);
t81 = pkin(8) + t53;
t31 = -t64 * t45 + t67 * t47;
t79 = t31 * t44;
t78 = t47 * t44;
t75 = t60 ^ 2 + t62 ^ 2;
t74 = 0.2e1 * t80;
t24 = t45 * t46;
t6 = t68 * t11 - t65 * t12;
t4 = t24 * pkin(9) + t6 + t85;
t1 = t67 * t4 - t64 * t5;
t41 = t81 * t60;
t42 = t81 * t62;
t27 = -t68 * t41 - t65 * t42;
t21 = pkin(4) * t77 + t34;
t72 = t16 * t62 + t17 * t60;
t71 = -t16 * t60 + t17 * t62;
t28 = -t65 * t41 + t68 * t42;
t70 = -t44 * t53 + t46 * t56;
t43 = t44 ^ 2;
t33 = t45 * t44;
t30 = t67 * t45 + t64 * t47;
t20 = t30 * t44;
t19 = -t45 * pkin(9) + t28;
t18 = -t47 * pkin(9) + t27;
t15 = t23 * pkin(5) + t21;
t14 = -t64 * t23 - t67 * t24;
t13 = t67 * t23 - t64 * t24;
t9 = t64 * t18 + t67 * t19;
t8 = t67 * t18 - t64 * t19;
t2 = t64 * t4 + t82;
t3 = [1, 0, 0, t66 ^ 2, t66 * t74, 0, 0, 0, pkin(1) * t74, -0.2e1 * pkin(1) * t66, 0.2e1 * t34 * t46 - 0.2e1 * t36 * t44, t36 ^ 2 + t57 ^ 2 + t90, 0.2e1 * t16 * t44 + 0.2e1 * t34 * t77, -0.2e1 * t17 * t44 + 0.2e1 * t34 * t76, -0.2e1 * t72 * t46, t16 ^ 2 + t17 ^ 2 + t90, t24 ^ 2, 0.2e1 * t24 * t23, -t24 * t87, t23 * t88, t43, 0.2e1 * t21 * t23 + 0.2e1 * t6 * t44, -0.2e1 * t21 * t24 - 0.2e1 * t7 * t44, t14 ^ 2, -0.2e1 * t14 * t13, t14 * t87, t13 * t88, t43, 0.2e1 * t1 * t44 + 0.2e1 * t15 * t13, 0.2e1 * t15 * t14 - 0.2e1 * t2 * t44; 0, 0, 0, 0, 0, t66, t80, 0, -t66 * pkin(7), -t73 (-t44 * t61 - t46 * t63) * pkin(2) (-t34 * t63 + t36 * t61) * pkin(2), -t34 * t62 + t60 * t70, t34 * t60 + t62 * t70, t71, t34 * t56 + t53 * t71, -t24 * t47, -t47 * t23 + t24 * t45, t78, -t33, 0, t21 * t45 + t49 * t23 + t27 * t44, t21 * t47 - t49 * t24 - t28 * t44, t14 * t31, -t31 * t13 - t14 * t30, t79, -t20, 0, t37 * t13 + t15 * t30 + t8 * t44, t37 * t14 + t15 * t31 - t9 * t44; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t61 ^ 2 + t63 ^ 2) * pkin(2) ^ 2, -0.2e1 * t56 * t62, 0.2e1 * t56 * t60, 0.2e1 * t75 * t53, t53 ^ 2 * t75 + t56 ^ 2, t47 ^ 2, -0.2e1 * t47 * t45, 0, 0, 0, t45 * t86, t47 * t86, t31 ^ 2, -0.2e1 * t31 * t30, 0, 0, 0, t30 * t89, t31 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t62 * t44, -t60 * t44, -t75 * t46, t72, 0, 0, 0, 0, 0, -t33, -t78, 0, 0, 0, 0, 0, -t20, -t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t75, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t76, 0, t34, 0, 0, 0, 0, 0, t23, -t24, 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t60, 0, t56, 0, 0, 0, 0, 0, t45, t47, 0, 0, 0, 0, 0, t30, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, t44, t6, -t7, 0, 0, t14, -t13, t44, t44 * t83 + t1, -t82 + (-t4 - t85) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t45, 0, t27, -t28, 0, 0, t31, -t30, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t47, 0, 0, 0, 0, 0, -t30, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t83, -0.2e1 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, t44, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, 0, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t83, -t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
