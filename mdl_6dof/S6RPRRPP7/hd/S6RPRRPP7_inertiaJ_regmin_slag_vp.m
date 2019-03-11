% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t38 = cos(qJ(4));
t56 = t38 * qJ(5);
t36 = sin(qJ(4));
t40 = pkin(4) + pkin(5);
t61 = t40 * t36;
t77 = t56 - t61;
t57 = t36 * qJ(5);
t76 = t38 * t40 + t57;
t75 = -0.2e1 * t36;
t74 = 0.2e1 * t38;
t39 = cos(qJ(3));
t73 = 0.2e1 * t39;
t72 = 0.2e1 * t40;
t71 = 2 * qJ(2);
t37 = sin(qJ(3));
t70 = pkin(8) * t37;
t69 = t36 * pkin(8);
t68 = t37 * pkin(4);
t24 = t36 * t37;
t67 = t36 * t38;
t25 = t36 * t39;
t42 = -pkin(1) - pkin(7);
t66 = t37 * t42;
t27 = t38 * t39;
t65 = t38 * t42;
t10 = pkin(3) + t76;
t64 = t39 * t10;
t49 = -t38 * pkin(4) - t57;
t17 = -pkin(3) + t49;
t63 = t39 * t17;
t62 = t39 * t37;
t28 = t39 * t42;
t16 = t37 * pkin(3) - t39 * pkin(8) + qJ(2);
t60 = -t38 * t16 + t36 * t66;
t7 = t36 * t16 + t37 * t65;
t32 = t36 ^ 2;
t34 = t38 ^ 2;
t59 = t32 + t34;
t33 = t37 ^ 2;
t35 = t39 ^ 2;
t58 = t33 + t35;
t55 = t38 * qJ(6);
t54 = -0.2e1 * t62;
t30 = t37 * qJ(5);
t4 = t30 + t7;
t53 = 0.2e1 * t30 + t7;
t14 = t59 * t37;
t52 = -pkin(3) * t39 - t70;
t5 = t60 - t68;
t51 = t5 * t36 + t4 * t38;
t50 = -t63 + t70;
t48 = pkin(4) * t36 - t56;
t18 = (pkin(8) - qJ(6)) * t36;
t31 = t38 * pkin(8);
t19 = t31 - t55;
t47 = t18 * t36 + t19 * t38;
t46 = -t39 * t55 + t60;
t44 = qJ(5) ^ 2;
t43 = 0.2e1 * qJ(5);
t26 = t38 * t37;
t23 = t37 * t56;
t22 = qJ(6) * t25;
t15 = t58 * t38;
t13 = t58 * t36;
t9 = t59 * t33 + t35;
t8 = t48 * t39 - t28;
t3 = t77 * t39 + t28;
t2 = t22 + t4;
t1 = -t40 * t37 + t46;
t6 = [1, 0, 0, -2 * pkin(1), t71, pkin(1) ^ 2 + qJ(2) ^ 2, t35, t54, 0, 0, 0, t37 * t71, t39 * t71, t34 * t35, -0.2e1 * t35 * t67, t62 * t74, t36 * t54, t33, -0.2e1 * t35 * t42 * t36 - 0.2e1 * t37 * t60, -0.2e1 * t35 * t65 - 0.2e1 * t7 * t37, 0.2e1 * t8 * t25 - 0.2e1 * t5 * t37 (-t36 * t4 + t38 * t5) * t73, -0.2e1 * t8 * t27 + 0.2e1 * t4 * t37, t4 ^ 2 + t5 ^ 2 + t8 ^ 2, -0.2e1 * t1 * t37 - 0.2e1 * t3 * t25, 0.2e1 * t2 * t37 + 0.2e1 * t3 * t27 (-t1 * t38 + t2 * t36) * t73, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, -t13, 0, t15, t51 * t37 - t8 * t39, -t13, t15, 0, t3 * t39 + (t1 * t36 + t2 * t38) * t37; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, 0, t28, -t66, t36 * t27 (-t32 + t34) * t39, t24, t26, 0, t38 * t28 + t52 * t36, -t36 * t28 + t52 * t38, -t50 * t36 - t8 * t38, t51, -t8 * t36 + t50 * t38, t51 * pkin(8) + t8 * t17, -t18 * t37 + t3 * t38 - t36 * t64, t19 * t37 + t3 * t36 + t38 * t64 (-t18 * t39 - t2) * t38 + (t19 * t39 - t1) * t36, t1 * t18 + t3 * t10 + t2 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, 0, 0, 0, 0, 0, t27, -t25, t27, t14, t25, pkin(8) * t14 - t63, t27, t25, -t14, t47 * t37 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t32, 0.2e1 * t67, 0, 0, 0, pkin(3) * t74, pkin(3) * t75, -0.2e1 * t17 * t38, 0.2e1 * t59 * pkin(8), t17 * t75, t59 * pkin(8) ^ 2 + t17 ^ 2, t10 * t74, 0.2e1 * t10 * t36, -0.2e1 * t47, t10 ^ 2 + t18 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, t37, -t60, -t7, -t60 + 0.2e1 * t68, t49 * t39, t53, -t5 * pkin(4) + t4 * qJ(5), t37 * t72 - t46, t22 + t53, t76 * t39, t2 * qJ(5) - t1 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26, -t24, 0, t26, -pkin(4) * t24 + t23, -t24, t26, 0, -t37 * t61 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t38, 0, -t69, -t31, -t69, -t48, t31, -t48 * pkin(8), -t18, t19, -t77, t19 * qJ(5) - t18 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, t43, pkin(4) ^ 2 + t44, t72, t43, 0, t40 ^ 2 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t27, 0, t5, -t37, 0, -t27, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t69, 0, 0, -t36, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), -1, 0, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t27, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t36, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
