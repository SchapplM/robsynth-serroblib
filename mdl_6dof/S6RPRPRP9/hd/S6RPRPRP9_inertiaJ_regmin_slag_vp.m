% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRP9
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
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t43 = cos(qJ(3));
t38 = sin(pkin(9));
t39 = cos(pkin(9));
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t74 = -t40 * t38 + t42 * t39;
t75 = t74 * t43;
t25 = t42 * t38 + t40 * t39;
t73 = -0.2e1 * t25;
t33 = -t39 * pkin(4) - pkin(3);
t72 = 0.2e1 * t33;
t41 = sin(qJ(3));
t71 = -0.2e1 * t41;
t70 = 2 * qJ(2);
t69 = t41 * pkin(5);
t68 = t43 * pkin(3);
t26 = t41 * pkin(3) - t43 * qJ(4) + qJ(2);
t44 = -pkin(1) - pkin(7);
t60 = t41 * t44;
t15 = t38 * t26 + t39 * t60;
t63 = t38 * t43;
t11 = -pkin(8) * t63 + t15;
t21 = t39 * t26;
t31 = t39 * t43;
t62 = t38 * t44;
t9 = -pkin(8) * t31 + t21 + (pkin(4) - t62) * t41;
t4 = t42 * t11 + t40 * t9;
t56 = pkin(8) + qJ(4);
t27 = t56 * t39;
t51 = t56 * t38;
t12 = t40 * t27 + t42 * t51;
t67 = t12 * t41;
t13 = t42 * t27 - t40 * t51;
t66 = t13 * t41;
t65 = t74 * t41;
t16 = t25 * t41;
t37 = t43 ^ 2;
t64 = t37 * t44;
t17 = t43 * t25;
t57 = t43 * t44;
t55 = t38 ^ 2 + t39 ^ 2;
t36 = t41 ^ 2;
t54 = -t36 - t37;
t53 = t41 * qJ(6);
t52 = t40 * t11 - t42 * t9;
t22 = pkin(4) * t63 - t57;
t50 = t55 * qJ(4);
t49 = -qJ(4) * t41 - t68;
t14 = -t38 * t60 + t21;
t48 = -t14 * t38 + t15 * t39;
t47 = -t16 * t41 - t43 * t17;
t46 = t41 * t65 + t43 * t75;
t8 = -pkin(5) * t74 - t25 * qJ(6) + t33;
t5 = t17 * pkin(5) - qJ(6) * t75 + t22;
t2 = t52 - t69;
t1 = t53 + t4;
t3 = [1, 0, 0, -2 * pkin(1), t70, pkin(1) ^ 2 + qJ(2) ^ 2, t37, t43 * t71, 0, 0, 0, t41 * t70, t43 * t70, 0.2e1 * t14 * t41 - 0.2e1 * t37 * t62, -0.2e1 * t15 * t41 - 0.2e1 * t39 * t64, 0.2e1 * (-t14 * t39 - t15 * t38) * t43, t37 * t44 ^ 2 + t14 ^ 2 + t15 ^ 2, t75 ^ 2, -0.2e1 * t75 * t17, 0.2e1 * t75 * t41, t17 * t71, t36, 0.2e1 * t22 * t17 - 0.2e1 * t41 * t52, 0.2e1 * t22 * t75 - 0.2e1 * t4 * t41, 0.2e1 * t5 * t17 - 0.2e1 * t2 * t41, -0.2e1 * t1 * t17 + 0.2e1 * t2 * t75, 0.2e1 * t1 * t41 - 0.2e1 * t5 * t75, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, t54 * t38, t54 * t39, 0, t48 * t41 + t64, 0, 0, 0, 0, 0, t47, -t46, t47, t16 * t75 - t17 * t65, t46, t1 * t65 + t2 * t16 - t5 * t43; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t36 + t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 ^ 2 + t65 ^ 2 + t37; 0, 0, 0, 0, 0, 0, 0, 0, t43, -t41, 0, t57, -t60, t49 * t38 + t39 * t57, -t38 * t57 + t49 * t39, t48, pkin(3) * t57 + t48 * qJ(4), t75 * t25, -t25 * t17 + t74 * t75, t16, t65, 0, t33 * t17 - t22 * t74 - t67, t22 * t25 + t33 * t75 - t66, t8 * t17 - t5 * t74 - t67, t1 * t74 + t12 * t75 - t13 * t17 + t2 * t25, -t5 * t25 - t75 * t8 + t66, t1 * t13 + t2 * t12 + t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t41, t31, -t63, t55 * t41, t41 * t50 + t68, 0, 0, 0, 0, 0, t75, -t17, t75, t16 * t25 + t65 * t74, t17, t16 * t12 + t13 * t65 - t43 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t39, -0.2e1 * pkin(3) * t38, 0.2e1 * t50, t55 * qJ(4) ^ 2 + pkin(3) ^ 2, t25 ^ 2, -t74 * t73, 0, 0, 0, -t74 * t72, t25 * t72, -0.2e1 * t8 * t74, 0.2e1 * t12 * t25 + 0.2e1 * t13 * t74, t8 * t73, t12 ^ 2 + t13 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t31, 0, -t57, 0, 0, 0, 0, 0, t17, t75, t17, 0, -t75, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, 0, -pkin(3), 0, 0, 0, 0, 0, -t74, t25, -t74, 0, -t25, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t17, t41, -t52, -t4, -t52 + 0.2e1 * t69, -pkin(5) * t75 - t17 * qJ(6), 0.2e1 * t53 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t65, -t16, 0, t65, -t16 * pkin(5) + qJ(6) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t74, 0, -t12, -t13, -t12, -pkin(5) * t25 + qJ(6) * t74, t13, -t12 * pkin(5) + t13 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t75, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
