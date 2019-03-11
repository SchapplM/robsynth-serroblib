% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t44 = sin(pkin(9));
t46 = cos(pkin(9));
t49 = sin(qJ(4));
t63 = cos(qJ(4));
t25 = t44 * t63 + t46 * t49;
t21 = t25 ^ 2;
t24 = t44 * t49 - t46 * t63;
t70 = t24 ^ 2;
t71 = -t21 - t70;
t43 = sin(pkin(10));
t45 = cos(pkin(10));
t48 = sin(qJ(6));
t50 = cos(qJ(6));
t26 = t43 * t50 + t45 * t48;
t62 = t24 * t26;
t69 = 0.2e1 * t62;
t68 = 0.2e1 * t25;
t36 = -pkin(5) * t45 - pkin(4);
t67 = 0.2e1 * t36;
t66 = 2 * qJ(2);
t65 = t24 * pkin(4);
t47 = -pkin(1) - qJ(3);
t64 = -pkin(7) + t47;
t33 = pkin(3) * t44 + qJ(2);
t14 = pkin(4) * t25 + qJ(5) * t24 + t33;
t28 = t64 * t44;
t29 = t64 * t46;
t17 = t28 * t63 + t29 * t49;
t6 = t14 * t43 + t17 * t45;
t23 = t43 * t48 - t45 * t50;
t20 = t23 * t25;
t13 = t24 * t23;
t61 = t24 * t43;
t60 = t24 * t45;
t59 = t26 * t25;
t58 = pkin(8) + qJ(5);
t57 = t43 ^ 2 + t45 ^ 2;
t32 = t44 ^ 2 + t46 ^ 2;
t5 = t14 * t45 - t17 * t43;
t56 = t57 * qJ(5);
t55 = t43 * t6 + t45 * t5;
t54 = -t43 * t5 + t45 * t6;
t53 = -qJ(5) * t25 + t65;
t16 = t28 * t49 - t29 * t63;
t52 = qJ(2) ^ 2;
t31 = t58 * t45;
t30 = t58 * t43;
t22 = t32 * t47;
t19 = -t30 * t48 + t31 * t50;
t18 = -t30 * t50 - t31 * t48;
t7 = -pkin(5) * t61 + t16;
t4 = pkin(8) * t61 + t6;
t3 = pkin(5) * t25 + pkin(8) * t60 + t5;
t2 = t3 * t48 + t4 * t50;
t1 = t3 * t50 - t4 * t48;
t8 = [1, 0, 0, -2 * pkin(1), t66, pkin(1) ^ 2 + t52, t44 * t66, t46 * t66, -0.2e1 * t22, t32 * t47 ^ 2 + t52, t70, t24 * t68, 0, 0, 0, t33 * t68, -0.2e1 * t33 * t24, -0.2e1 * t16 * t61 + 0.2e1 * t25 * t5, -0.2e1 * t16 * t60 - 0.2e1 * t25 * t6, 0.2e1 * t55 * t24, t16 ^ 2 + t5 ^ 2 + t6 ^ 2, t13 ^ 2, t13 * t69, t13 * t68, t25 * t69, t21, 0.2e1 * t1 * t25 - 0.2e1 * t62 * t7, 0.2e1 * t13 * t7 - 0.2e1 * t2 * t25; 0, 0, 0, 1, 0, -pkin(1), 0, 0, -t32, t22, 0, 0, 0, 0, 0, 0, 0, t71 * t43, t71 * t45, 0, t16 * t24 + t25 * t54, 0, 0, 0, 0, 0, -t24 * t62 - t25 * t59, t13 * t24 + t20 * t25; 0, 0, 0, 0, 0, 1, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t57 + t70, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t44, t46, 0, qJ(2), 0, 0, 0, 0, 0, t25, -t24, t45 * t25, -t43 * t25, t57 * t24, t55, 0, 0, 0, 0, 0, -t20, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, 0, -t16, -t17, -t16 * t45 + t43 * t53, t16 * t43 + t45 * t53, t54, -t16 * pkin(4) + qJ(5) * t54, t13 * t26, -t13 * t23 + t26 * t62, t59, -t20, 0, t18 * t25 + t23 * t7 - t36 * t62, t13 * t36 - t19 * t25 + t26 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, -t60, t61, t57 * t25, t25 * t56 - t65, 0, 0, 0, 0, 0, t13, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t45, -0.2e1 * pkin(4) * t43, 0.2e1 * t56, qJ(5) ^ 2 * t57 + pkin(4) ^ 2, t26 ^ 2, -0.2e1 * t26 * t23, 0, 0, 0, t23 * t67, t26 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t60, 0, t16, 0, 0, 0, 0, 0, -t62, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t43, 0, -pkin(4), 0, 0, 0, 0, 0, t23, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t62, t25, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t23, 0, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
