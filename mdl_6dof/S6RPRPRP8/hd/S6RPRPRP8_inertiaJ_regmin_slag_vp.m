% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t39 = cos(qJ(5));
t36 = sin(pkin(9));
t37 = cos(pkin(9));
t40 = cos(qJ(3));
t63 = sin(qJ(3));
t21 = -t36 * t63 + t37 * t40;
t19 = t21 ^ 2;
t22 = -t36 * t40 - t37 * t63;
t20 = t22 ^ 2;
t68 = t20 + t19;
t71 = t68 * t39;
t70 = -0.2e1 * t21;
t69 = (t21 * t37 - t22 * t36) * pkin(3);
t38 = sin(qJ(5));
t67 = t68 * t38;
t66 = -0.2e1 * t39;
t65 = 2 * qJ(2);
t41 = -pkin(1) - pkin(7);
t53 = t63 * t41;
t25 = -t63 * qJ(4) + t53;
t51 = (-qJ(4) + t41) * t40;
t12 = t37 * t25 + t36 * t51;
t32 = t63 * pkin(3) + qJ(2);
t8 = -t22 * pkin(4) - t21 * pkin(8) + t32;
t4 = t39 * t12 + t38 * t8;
t64 = pkin(5) * t22;
t31 = -t37 * pkin(3) - pkin(4);
t48 = t39 * pkin(5) + t38 * qJ(6);
t17 = t31 - t48;
t62 = t21 * t17;
t30 = t36 * pkin(3) + pkin(8);
t61 = t22 * t30;
t60 = t38 * t21;
t14 = t38 * t22;
t59 = t38 * t30;
t58 = t38 * t39;
t15 = t39 * t21;
t16 = t39 * t22;
t57 = t39 * t30;
t34 = t38 ^ 2;
t35 = t39 ^ 2;
t56 = t34 + t35;
t55 = t22 * qJ(6);
t3 = -t38 * t12 + t39 * t8;
t52 = t56 * t22;
t10 = t36 * t25 - t37 * t51;
t1 = -t55 + t4;
t2 = -t3 + t64;
t50 = t1 * t39 + t2 * t38;
t49 = t1 * t38 - t2 * t39;
t47 = pkin(5) * t38 - t39 * qJ(6);
t46 = t10 * t21 + t12 * t22;
t45 = t61 + t62;
t44 = t21 * t31 + t61;
t5 = t47 * t21 + t10;
t6 = [1, 0, 0, -2 * pkin(1), t65, pkin(1) ^ 2 + qJ(2) ^ 2, t40 ^ 2, -0.2e1 * t40 * t63, 0, 0, 0, t63 * t65, t40 * t65, 0.2e1 * t46, t10 ^ 2 + t12 ^ 2 + t32 ^ 2, t35 * t19, -0.2e1 * t19 * t58, t16 * t70, 0.2e1 * t21 * t14, t20, 0.2e1 * t10 * t60 - 0.2e1 * t3 * t22, 0.2e1 * t10 * t15 + 0.2e1 * t4 * t22, 0.2e1 * t2 * t22 + 0.2e1 * t5 * t60, t49 * t70, -0.2e1 * t1 * t22 - 0.2e1 * t5 * t15, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t68, -t46, 0, 0, 0, 0, 0, -t67, -t71, -t67, 0, t71, -t5 * t21 - t50 * t22; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t20 + t19; 0, 0, 0, 0, 0, 0, 0, 0, t40, -t63, 0, t40 * t41, -t53, -t69 (-t10 * t37 + t12 * t36) * pkin(3), t38 * t15 (-t34 + t35) * t21, -t14, -t16, 0, -t10 * t39 + t44 * t38, t10 * t38 + t44 * t39, t45 * t38 - t5 * t39, t50, -t5 * t38 - t45 * t39, t5 * t17 + t50 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t63, 0, t69, 0, 0, 0, 0, 0, t15, -t60, t15, -t52, t60, -t30 * t52 - t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t36 ^ 2 + t37 ^ 2) * pkin(3) ^ 2, t34, 0.2e1 * t58, 0, 0, 0, t31 * t66, 0.2e1 * t31 * t38, t17 * t66, 0.2e1 * t56 * t30, -0.2e1 * t17 * t38, t56 * t30 ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, -t16, t14, -t16, -t56 * t21, -t14, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t60, -t22, t3, -t4, t3 - 0.2e1 * t64, -t48 * t21, -0.2e1 * t55 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, t14, 0, -t16, t47 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t39, 0, -t59, -t57, -t59, -t47, t57, -t47 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t38, t39, 0, t38, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t15, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
