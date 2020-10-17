% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:37:38
% EndTime: 2019-05-05 06:37:39
% DurationCPUTime: 0.57s
% Computational Cost: add. (812->118), mult. (1747->215), div. (0->0), fcn. (1988->10), ass. (0->65)
t55 = sin(qJ(3));
t86 = -0.2e1 * t55;
t58 = cos(qJ(3));
t85 = 0.2e1 * t58;
t57 = cos(qJ(4));
t84 = pkin(3) * t57;
t54 = sin(qJ(4));
t83 = pkin(8) * t54;
t52 = sin(pkin(6));
t82 = t52 * sin(qJ(2));
t81 = t52 * cos(qJ(2));
t80 = t54 * t55;
t79 = t54 * t57;
t78 = t54 * t58;
t77 = t57 * t55;
t76 = t57 * t58;
t75 = -qJ(5) - pkin(9);
t38 = -t58 * pkin(3) - t55 * pkin(9) - pkin(2);
t35 = t57 * t38;
t74 = qJ(5) * t55;
t19 = -t57 * t74 + t35 + (-pkin(4) - t83) * t58;
t70 = pkin(8) * t76;
t22 = t70 + (t38 - t74) * t54;
t51 = sin(pkin(11));
t53 = cos(pkin(11));
t7 = t51 * t19 + t53 * t22;
t47 = t55 * pkin(8);
t37 = pkin(4) * t80 + t47;
t73 = cos(pkin(6));
t72 = t55 * t85;
t39 = t75 * t57;
t65 = t75 * t54;
t24 = -t51 * t39 - t53 * t65;
t26 = -t53 * t39 + t51 * t65;
t71 = t24 ^ 2 + t26 ^ 2;
t46 = -t57 * pkin(4) - pkin(3);
t62 = t73 * t55 + t58 * t82;
t21 = -t54 * t81 + t62 * t57;
t61 = t62 * t54 + t57 * t81;
t11 = t53 * t21 - t51 * t61;
t9 = t51 * t21 + t53 * t61;
t69 = t11 * t26 + t9 * t24;
t34 = t51 * t57 + t53 * t54;
t29 = t34 * t55;
t30 = -t51 * t80 + t53 * t77;
t68 = -t11 * t29 + t9 * t30;
t33 = t51 * t54 - t53 * t57;
t67 = -t11 * t33 + t9 * t34;
t32 = t55 * t82 - t73 * t58;
t66 = t11 ^ 2 + t32 ^ 2 + t9 ^ 2;
t64 = t24 * t30 - t26 * t29;
t6 = t53 * t19 - t51 * t22;
t63 = 0.2e1 * t24 * t34 - 0.2e1 * t26 * t33;
t50 = t57 ^ 2;
t49 = t55 ^ 2;
t48 = t54 ^ 2;
t44 = t53 * pkin(4) + pkin(5);
t42 = t51 * pkin(4) + qJ(6);
t28 = t54 * t38 + t70;
t27 = -pkin(8) * t78 + t35;
t15 = t33 * pkin(5) - t34 * qJ(6) + t46;
t12 = t29 * pkin(5) - t30 * qJ(6) + t37;
t5 = t58 * pkin(5) - t6;
t4 = -t58 * qJ(6) + t7;
t1 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, t66; 0, 0, t81, -t82, 0, 0, 0, 0, 0, t58 * t81, -t55 * t81, 0, 0, 0, 0, 0, t32 * t80 + t61 * t58, t21 * t58 + t32 * t77, t68, t11 * t7 + t32 * t37 - t9 * t6, t32 * t29 + t9 * t58, t68, -t11 * t58 - t32 * t30, t11 * t4 + t32 * t12 + t9 * t5; 0, 1, 0, 0, t49, t72, 0, 0, 0, pkin(2) * t85, pkin(2) * t86, t50 * t49, -0.2e1 * t49 * t79, t76 * t86, t54 * t72, t58 ^ 2, -0.2e1 * t27 * t58 + 0.2e1 * t49 * t83, 0.2e1 * t49 * pkin(8) * t57 + 0.2e1 * t28 * t58, -0.2e1 * t7 * t29 - 0.2e1 * t6 * t30, t37 ^ 2 + t6 ^ 2 + t7 ^ 2, 0.2e1 * t12 * t29 + 0.2e1 * t5 * t58, -0.2e1 * t4 * t29 + 0.2e1 * t5 * t30, -0.2e1 * t12 * t30 - 0.2e1 * t4 * t58, t12 ^ 2 + t4 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t62, 0, 0, 0, 0, 0, -t32 * t57, t32 * t54, t67, t32 * t46 + t69, t32 * t33, t67, -t32 * t34, t32 * t15 + t69; 0, 0, 0, 0, 0, 0, t55, t58, 0, -t47, -t58 * pkin(8), t54 * t77 (-t48 + t50) * t55, -t78, -t76, 0, -pkin(8) * t77 + (-pkin(3) * t55 + pkin(9) * t58) * t54, pkin(9) * t76 + (t83 - t84) * t55, -t7 * t33 - t6 * t34 + t64, -t6 * t24 + t7 * t26 + t37 * t46, t12 * t33 + t15 * t29 + t24 * t58, -t4 * t33 + t5 * t34 + t64, -t12 * t34 - t15 * t30 - t26 * t58, t12 * t15 + t5 * t24 + t4 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t48, 0.2e1 * t79, 0, 0, 0, 0.2e1 * t84, -0.2e1 * pkin(3) * t54, t63, t46 ^ 2 + t71, 0.2e1 * t15 * t33, t63, -0.2e1 * t15 * t34, t15 ^ 2 + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t21, 0 (t11 * t51 - t53 * t9) * pkin(4), -t9, 0, t11, t11 * t42 - t9 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t80, -t58, t27, -t28 (-t29 * t51 - t30 * t53) * pkin(4) (t51 * t7 + t53 * t6) * pkin(4) (-pkin(5) - t44) * t58 + t6, -t42 * t29 - t44 * t30 (-qJ(6) - t42) * t58 + t7, t4 * t42 - t5 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t57, 0, -t54 * pkin(9), -t57 * pkin(9) (-t33 * t51 - t34 * t53) * pkin(4) (-t24 * t53 + t26 * t51) * pkin(4), -t24, -t42 * t33 - t44 * t34, t26, -t24 * t44 + t26 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t51 ^ 2 + t53 ^ 2) * pkin(4) ^ 2, 0.2e1 * t44, 0, 0.2e1 * t42, t42 ^ 2 + t44 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t29, 0, -t30, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t33, 0, -t34, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t30, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t1;
