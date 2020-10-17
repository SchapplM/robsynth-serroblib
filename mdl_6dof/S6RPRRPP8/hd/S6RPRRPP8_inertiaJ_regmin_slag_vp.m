% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPP8
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:52:54
% EndTime: 2019-05-05 21:52:56
% DurationCPUTime: 0.66s
% Computational Cost: add. (423->100), mult. (732->173), div. (0->0), fcn. (684->4), ass. (0->69)
t40 = cos(qJ(4));
t29 = qJ(5) * t40;
t38 = sin(qJ(4));
t77 = -pkin(4) * t38 + t29;
t36 = pkin(4) + qJ(6);
t56 = t38 * qJ(5);
t76 = -t36 * t40 - t56;
t75 = 0.2e1 * t36;
t74 = -0.2e1 * t38;
t73 = 0.2e1 * t40;
t41 = cos(qJ(3));
t72 = 0.2e1 * t41;
t71 = 2 * qJ(2);
t39 = sin(qJ(3));
t69 = pkin(8) * t39;
t68 = t39 * pkin(4);
t24 = t38 * t39;
t66 = t38 * t40;
t25 = t38 * t41;
t42 = -pkin(1) - pkin(7);
t65 = t39 * t42;
t27 = t40 * t41;
t64 = t40 * t42;
t10 = -pkin(3) + t76;
t63 = t41 * t10;
t49 = -t40 * pkin(4) - t56;
t17 = -pkin(3) + t49;
t62 = t41 * t17;
t61 = t41 * t39;
t60 = t41 * t42;
t16 = t39 * pkin(3) - t41 * pkin(8) + qJ(2);
t59 = -t40 * t16 + t38 * t65;
t7 = t38 * t16 + t39 * t64;
t32 = t38 ^ 2;
t34 = t40 ^ 2;
t58 = t32 + t34;
t33 = t39 ^ 2;
t35 = t41 ^ 2;
t57 = t33 + t35;
t55 = t39 * qJ(5);
t54 = -0.2e1 * t61;
t14 = t58 * t39;
t53 = -t42 - t29;
t52 = -pkin(3) * t41 - t69;
t4 = -t55 - t7;
t5 = t59 - t68;
t51 = t5 * t38 - t4 * t40;
t50 = -t62 + t69;
t30 = t38 * pkin(8);
t18 = t38 * pkin(5) + t30;
t31 = t40 * pkin(8);
t19 = t40 * pkin(5) + t31;
t48 = t18 * t38 + t19 * t40;
t47 = -pkin(5) * t25 + t7;
t46 = -pkin(5) * t27 - t59;
t44 = qJ(5) ^ 2;
t43 = 0.2e1 * qJ(5);
t28 = 0.2e1 * t55;
t26 = t40 * t39;
t23 = pkin(4) * t25;
t22 = t40 * t55;
t15 = t57 * t40;
t13 = t57 * t38;
t9 = t58 * t33 + t35;
t8 = t53 * t41 + t23;
t3 = t23 + (qJ(6) * t38 + t53) * t41;
t2 = t47 + t55;
t1 = -t36 * t39 - t46;
t6 = [1, 0, 0, -2 * pkin(1), t71, pkin(1) ^ 2 + qJ(2) ^ 2, t35, t54, 0, 0, 0, t39 * t71, t41 * t71, t34 * t35, -0.2e1 * t35 * t66, t61 * t73, t38 * t54, t33, -0.2e1 * t35 * t42 * t38 - 0.2e1 * t39 * t59, -0.2e1 * t35 * t64 - 0.2e1 * t7 * t39 (t38 * t4 + t40 * t5) * t72, -0.2e1 * t8 * t25 + 0.2e1 * t5 * t39, -0.2e1 * t8 * t27 - 0.2e1 * t4 * t39, t4 ^ 2 + t5 ^ 2 + t8 ^ 2 (t1 * t40 - t2 * t38) * t72, 0.2e1 * t2 * t39 - 0.2e1 * t3 * t27, -0.2e1 * t1 * t39 + 0.2e1 * t3 * t25, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, t13, t15, t51 * t39 - t8 * t41, 0, t15, -t13, -t3 * t41 + (t1 * t38 + t2 * t40) * t39; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, t41, -t39, 0, t60, -t65, t38 * t27 (-t32 + t34) * t41, t24, t26, 0, t52 * t38 + t40 * t60, -t38 * t60 + t52 * t40, t51, t50 * t38 + t8 * t40, -t8 * t38 + t50 * t40, t51 * pkin(8) + t8 * t17 (t18 * t41 + t2) * t40 + (-t19 * t41 + t1) * t38, t19 * t39 - t3 * t38 - t40 * t63, -t18 * t39 - t3 * t40 + t38 * t63, t1 * t18 + t3 * t10 + t2 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t39, 0, 0, 0, 0, 0, t27, -t25, t14, -t27, t25, pkin(8) * t14 - t62, t14, t25, t27, t48 * t39 - t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t32, 0.2e1 * t66, 0, 0, 0, pkin(3) * t73, pkin(3) * t74, 0.2e1 * t58 * pkin(8), t17 * t73, t17 * t74, t58 * pkin(8) ^ 2 + t17 ^ 2, 0.2e1 * t48, t10 * t74, -0.2e1 * t10 * t40, t10 ^ 2 + t18 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, t39, -t59, -t7, t49 * t41, t59 - 0.2e1 * t68, t28 + t7, -t5 * pkin(4) - t4 * qJ(5), t76 * t41, t28 + t47, t39 * t75 + t46, t2 * qJ(5) - t1 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t26, 0, t24, t26, -pkin(4) * t24 + t22, 0, t26, -t24, -t36 * t24 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t40, 0, -t30, -t31, t77, t30, t31, t77 * pkin(8), -t36 * t38 + t29, t19, -t18, t19 * qJ(5) - t18 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t43, pkin(4) ^ 2 + t44, 0, t43, t75, t36 ^ 2 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t39, 0, t5, t27, 0, -t39, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, t30, t38, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, -1, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t39, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
