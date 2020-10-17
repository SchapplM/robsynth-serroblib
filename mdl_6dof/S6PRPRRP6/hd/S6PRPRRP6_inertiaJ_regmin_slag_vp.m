% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:02:08
% EndTime: 2019-05-05 00:02:09
% DurationCPUTime: 0.50s
% Computational Cost: add. (295->88), mult. (627->155), div. (0->0), fcn. (689->8), ass. (0->63)
t68 = 2 * pkin(5);
t30 = sin(qJ(5));
t67 = -0.2e1 * t30;
t33 = cos(qJ(5));
t66 = 0.2e1 * t33;
t65 = 2 * qJ(3);
t31 = sin(qJ(4));
t64 = pkin(9) * t31;
t63 = t30 * pkin(9);
t62 = t33 * pkin(9);
t29 = cos(pkin(6));
t34 = cos(qJ(4));
t28 = sin(pkin(6));
t35 = cos(qJ(2));
t58 = t28 * t35;
t10 = t29 * t31 + t34 * t58;
t61 = t10 * t30;
t60 = t10 * t33;
t32 = sin(qJ(2));
t59 = t28 * t32;
t57 = t30 * t33;
t56 = t30 * t34;
t36 = -pkin(2) - pkin(8);
t55 = t30 * t36;
t54 = t31 * t36;
t22 = t33 * t34;
t53 = t33 * t36;
t39 = -t33 * pkin(5) - t30 * qJ(6);
t17 = -pkin(4) + t39;
t52 = t34 * t17;
t51 = t34 * t31;
t50 = t34 * t36;
t16 = t31 * pkin(4) - t34 * pkin(9) + qJ(3);
t7 = t30 * t16 + t31 * t53;
t24 = t30 ^ 2;
t26 = t33 ^ 2;
t49 = t24 + t26;
t25 = t31 ^ 2;
t27 = t34 ^ 2;
t48 = t25 + t27;
t47 = t31 * qJ(6);
t46 = -0.2e1 * t51;
t11 = t29 * t34 - t31 * t58;
t2 = t11 * t30 - t33 * t59;
t45 = t10 * t56 - t2 * t31;
t44 = t49 * t31;
t43 = -pkin(4) * t34 - t64;
t3 = t11 * t33 + t30 * t59;
t42 = t2 * t30 + t3 * t33;
t4 = t47 + t7;
t13 = t33 * t16;
t5 = -t13 + (-pkin(5) + t55) * t31;
t41 = t5 * t30 + t4 * t33;
t40 = -t52 + t64;
t38 = -pkin(5) * t30 + t33 * qJ(6);
t21 = t33 * t31;
t20 = t30 * t31;
t15 = t48 * t33;
t14 = t48 * t30;
t9 = (-t36 - t38) * t34;
t6 = -t30 * t54 + t13;
t1 = -t10 * t22 + t3 * t31;
t8 = [1, 0, 0, 0, 0, 0, t29 ^ 2 + (t32 ^ 2 + t35 ^ 2) * t28 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, t58, -t59, -t58, t59 (pkin(2) * t35 + qJ(3) * t32) * t28, 0, 0, 0, 0, 0, t31 * t59, t34 * t59, 0, 0, 0, 0, 0, t45, -t1, t45 (t2 * t33 - t3 * t30) * t34, t1, t10 * t9 + t2 * t5 + t3 * t4; 0, 1, 0, 0, -0.2e1 * pkin(2), t65, pkin(2) ^ 2 + (qJ(3) ^ 2) t27, t46, 0, 0, 0, t31 * t65, t34 * t65, t26 * t27, -0.2e1 * t27 * t57, t51 * t66, t30 * t46, t25, -0.2e1 * t27 * t55 + 0.2e1 * t6 * t31, -0.2e1 * t27 * t53 - 0.2e1 * t7 * t31, -0.2e1 * t5 * t31 + 0.2e1 * t9 * t56, 0.2e1 * (-t30 * t4 + t33 * t5) * t34, -0.2e1 * t9 * t22 + 0.2e1 * t4 * t31, t4 ^ 2 + t5 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, -t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t34 + t42 * t31; 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, -t14, 0, t15, t41 * t31 - t9 * t34; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t25 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t60, t61, -t60, t42, -t61, t42 * pkin(9) + t10 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t31, 0, t50, -t54, t30 * t22 (-t24 + t26) * t34, t20, t21, 0, t43 * t30 + t33 * t50, -t30 * t50 + t43 * t33, -t40 * t30 - t9 * t33, t41, -t9 * t30 + t40 * t33, t41 * pkin(9) + t9 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t31, 0, 0, 0, 0, 0, t22, -t56, t22, t44, t56, pkin(9) * t44 - t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t24, 0.2e1 * t57, 0, 0, 0, pkin(4) * t66, pkin(4) * t67, -0.2e1 * t17 * t33, 0.2e1 * t49 * pkin(9), t17 * t67, t49 * pkin(9) ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, -t3, -t2, 0, t3, -t2 * pkin(5) + t3 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t56, t31, t6, -t7, t13 + (t68 - t55) * t31, t39 * t34, 0.2e1 * t47 + t7, -t5 * pkin(5) + t4 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, -t20, 0, t21, t38 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t33, 0, -t63, -t62, -t63, t38, t62, t38 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t68, 0, 0.2e1 * qJ(6) (pkin(5) ^ 2) + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t22, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t8;
