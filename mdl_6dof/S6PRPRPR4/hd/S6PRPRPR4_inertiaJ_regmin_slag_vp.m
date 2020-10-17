% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRPR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:47:51
% EndTime: 2019-05-04 22:47:52
% DurationCPUTime: 0.56s
% Computational Cost: add. (613->97), mult. (1356->183), div. (0->0), fcn. (1698->12), ass. (0->63)
t44 = sin(pkin(11));
t47 = cos(pkin(11));
t50 = sin(qJ(4));
t71 = cos(qJ(4));
t28 = t50 * t44 - t71 * t47;
t74 = -0.2e1 * t28;
t46 = cos(pkin(12));
t37 = -t46 * pkin(5) - pkin(4);
t73 = 0.2e1 * t37;
t38 = -t47 * pkin(3) - pkin(2);
t72 = 0.2e1 * t38;
t43 = sin(pkin(12));
t49 = sin(qJ(6));
t52 = cos(qJ(6));
t29 = t52 * t43 + t49 * t46;
t70 = t29 * t28;
t30 = t71 * t44 + t50 * t47;
t69 = t43 * t30;
t45 = sin(pkin(6));
t68 = t45 * sin(qJ(2));
t53 = cos(qJ(2));
t67 = t45 * t53;
t66 = t46 * t30;
t65 = pkin(8) + qJ(3);
t64 = pkin(9) + qJ(5);
t18 = t28 * pkin(4) - t30 * qJ(5) + t38;
t32 = t65 * t44;
t34 = t65 * t47;
t23 = -t50 * t32 + t71 * t34;
t8 = t43 * t18 + t46 * t23;
t63 = t43 ^ 2 + t46 ^ 2;
t62 = t44 ^ 2 + t47 ^ 2;
t7 = t46 * t18 - t43 * t23;
t61 = t8 * t43 + t7 * t46;
t60 = -t7 * t43 + t8 * t46;
t48 = cos(pkin(6));
t25 = -t44 * t68 + t48 * t47;
t26 = t48 * t44 + t47 * t68;
t13 = t50 * t25 + t71 * t26;
t10 = t46 * t13 - t43 * t67;
t9 = -t43 * t13 - t46 * t67;
t59 = t10 * t46 - t9 * t43;
t58 = t10 * t43 + t9 * t46;
t57 = -pkin(4) * t30 - qJ(5) * t28;
t56 = -t25 * t44 + t26 * t47;
t27 = t49 * t43 - t52 * t46;
t21 = t71 * t32 + t50 * t34;
t33 = t64 * t46;
t31 = t64 * t43;
t24 = t27 * t28;
t22 = -t49 * t31 + t52 * t33;
t20 = -t52 * t31 - t49 * t33;
t15 = t27 * t30;
t14 = t29 * t30;
t12 = -t71 * t25 + t50 * t26;
t11 = pkin(5) * t69 + t21;
t6 = -pkin(9) * t69 + t8;
t5 = t28 * pkin(5) - pkin(9) * t66 + t7;
t4 = t52 * t10 + t49 * t9;
t3 = -t49 * t10 + t52 * t9;
t2 = t49 * t5 + t52 * t6;
t1 = -t49 * t6 + t52 * t5;
t16 = [1, 0, 0, 0, 0, 0, 0, t45 ^ 2 * t53 ^ 2 + t25 ^ 2 + t26 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t12 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t67, -t68, t47 * t67, -t44 * t67, t56, pkin(2) * t67 + t56 * qJ(3), 0, 0, 0, 0, 0, -t28 * t67, -t30 * t67, t12 * t69 + t9 * t28, -t10 * t28 + t12 * t66, -t58 * t30, t10 * t8 + t12 * t21 + t9 * t7, 0, 0, 0, 0, 0, t12 * t14 + t3 * t28, -t12 * t15 - t4 * t28; 0, 1, 0, 0, 0.2e1 * pkin(2) * t47, -0.2e1 * pkin(2) * t44, 0.2e1 * t62 * qJ(3), t62 * qJ(3) ^ 2 + pkin(2) ^ 2, t30 ^ 2, t30 * t74, 0, 0, 0, t28 * t72, t30 * t72, 0.2e1 * t21 * t69 + 0.2e1 * t7 * t28, 0.2e1 * t21 * t66 - 0.2e1 * t8 * t28, -0.2e1 * t61 * t30, t21 ^ 2 + t7 ^ 2 + t8 ^ 2, t15 ^ 2, 0.2e1 * t15 * t14, t15 * t74, t14 * t74, t28 ^ 2, 0.2e1 * t1 * t28 + 0.2e1 * t11 * t14, -0.2e1 * t11 * t15 - 0.2e1 * t2 * t28; 0, 0, 0, 0, 0, 0, 0, -t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t47, t44, 0, -pkin(2), 0, 0, 0, 0, 0, t28, t30, t46 * t28, -t43 * t28, -t63 * t30, t61, 0, 0, 0, 0, 0, -t24, -t70; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, -t12 * t46, t12 * t43, t59, -t12 * pkin(4) + t59 * qJ(5), 0, 0, 0, 0, 0, t12 * t27, t12 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t28, 0, -t21, -t23, -t21 * t46 + t57 * t43, t21 * t43 + t57 * t46, t60, -t21 * pkin(4) + t60 * qJ(5), -t15 * t29, -t29 * t14 + t15 * t27, t70, -t24, 0, t11 * t27 + t37 * t14 + t20 * t28, t11 * t29 - t37 * t15 - t22 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t46, -0.2e1 * pkin(4) * t43, 0.2e1 * t63 * qJ(5), t63 * qJ(5) ^ 2 + pkin(4) ^ 2, t29 ^ 2, -0.2e1 * t29 * t27, 0, 0, 0, t27 * t73, t29 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t66, 0, t21, 0, 0, 0, 0, 0, t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t43, 0, -pkin(4), 0, 0, 0, 0, 0, t27, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t14, t28, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t27, 0, t20, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t16;
