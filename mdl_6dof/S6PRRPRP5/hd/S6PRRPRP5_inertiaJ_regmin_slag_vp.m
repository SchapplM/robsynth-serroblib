% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:10:55
% EndTime: 2019-05-05 04:10:57
% DurationCPUTime: 0.60s
% Computational Cost: add. (350->92), mult. (707->159), div. (0->0), fcn. (773->8), ass. (0->65)
t42 = cos(pkin(6));
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t41 = sin(pkin(6));
t71 = t41 * sin(qJ(2));
t18 = t42 * t44 + t47 * t71;
t77 = t18 ^ 2;
t76 = -0.2e1 * t44;
t75 = 0.2e1 * t47;
t74 = 2 * qJ(4);
t49 = -pkin(3) - pkin(9);
t73 = t44 * pkin(5);
t46 = cos(qJ(5));
t72 = t18 * t46;
t48 = cos(qJ(2));
t70 = t41 * t48;
t43 = sin(qJ(5));
t69 = t43 * t44;
t68 = t43 * t47;
t67 = t43 * t49;
t66 = t44 * t47;
t65 = t44 * t49;
t64 = t46 * t43;
t63 = t46 * t47;
t32 = t46 * t49;
t55 = -t44 * qJ(4) - pkin(2);
t20 = t49 * t47 + t55;
t33 = t44 * pkin(8);
t27 = t44 * pkin(4) + t33;
t7 = t46 * t20 + t43 * t27;
t34 = t47 * pkin(8);
t28 = t47 * pkin(4) + t34;
t36 = t43 ^ 2;
t38 = t46 ^ 2;
t30 = t36 + t38;
t37 = t44 ^ 2;
t39 = t47 ^ 2;
t62 = t37 + t39;
t61 = qJ(4) * t47;
t60 = t44 * qJ(6);
t59 = -0.2e1 * t66;
t58 = t44 * t70;
t57 = t47 * t70;
t17 = -t42 * t47 + t44 * t71;
t8 = t17 * t46 + t43 * t70;
t56 = t18 * t63 + t8 * t44;
t54 = t43 * t20 - t46 * t27;
t4 = t60 + t7;
t5 = t54 - t73;
t1 = t4 * t43 - t5 * t46;
t10 = -t17 * t43 + t46 * t70;
t2 = -t10 * t43 + t8 * t46;
t53 = -pkin(3) * t44 + t61;
t26 = pkin(5) * t46 + t43 * qJ(6);
t52 = t43 * pkin(5) - t46 * qJ(6);
t51 = t17 * t44 + t18 * t47;
t31 = t46 * t44;
t29 = t44 * t32;
t25 = -t47 * pkin(3) + t55;
t24 = qJ(4) + t52;
t23 = t30 * t49;
t15 = t18 * t43;
t12 = t26 * t47 + t28;
t3 = -t10 * t44 + t18 * t68;
t6 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 ^ 2 * t48 ^ 2 + t17 ^ 2 + t77, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t8 ^ 2 + t77; 0, 0, t70, -t71, 0, 0, 0, 0, 0, t57, -t58, t51, -t57, t58, t51 * pkin(8) - t25 * t70, 0, 0, 0, 0, 0, t56, -t3, t56 (t10 * t46 + t43 * t8) * t47, t3, -t10 * t4 + t18 * t12 - t8 * t5; 0, 1, 0, 0, t37, 0.2e1 * t66, 0, 0, 0, pkin(2) * t75, pkin(2) * t76, 0.2e1 * t62 * pkin(8), t25 * t75, t25 * t76, t62 * pkin(8) ^ 2 + t25 ^ 2, t36 * t39, 0.2e1 * t39 * t64, t43 * t59, t46 * t59, t37, 0.2e1 * t28 * t63 - 0.2e1 * t44 * t54, -0.2e1 * t28 * t68 - 0.2e1 * t7 * t44, 0.2e1 * t12 * t63 - 0.2e1 * t5 * t44 (-t4 * t46 - t43 * t5) * t75, 0.2e1 * t12 * t68 + 0.2e1 * t4 * t44, t12 ^ 2 + t4 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t18, 0, t17, t18, -t17 * pkin(3) + t18 * qJ(4), 0, 0, 0, 0, 0, t15, t72, t15, -t2, -t72, t18 * t24 + t2 * t49; 0, 0, 0, 0, 0, 0, t44, t47, 0, -t33, -t34, t53, t33, t34, t53 * pkin(8), -t43 * t63 (t36 - t38) * t47, t31, -t69, 0, t28 * t43 + t46 * t61 + t29, t28 * t46 + (-t61 - t65) * t43, t12 * t43 + t24 * t63 + t29, -t1, -t12 * t46 + (t24 * t47 + t65) * t43, t1 * t49 + t12 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t74, pkin(3) ^ 2 + (qJ(4) ^ 2) t38, -0.2e1 * t64, 0, 0, 0, t43 * t74, t46 * t74, 0.2e1 * t24 * t43, -0.2e1 * t23, -0.2e1 * t24 * t46, t30 * t49 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, t33, 0, 0, 0, 0, 0, t31, -t69, t31, 0, t69, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t10, t8, 0, -t10, t8 * pkin(5) - t10 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t63, t44, -t54, -t7, -t54 + 0.2e1 * t73, t52 * t47, 0.2e1 * t60 + t7, -t5 * pkin(5) + t4 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t43, 0, t32, -t67, t32, -t26, t67, t26 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t43, t46, 0, t43, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t68, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
