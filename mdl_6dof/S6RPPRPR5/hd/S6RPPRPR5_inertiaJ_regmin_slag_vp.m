% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:21:54
% EndTime: 2019-05-05 14:21:56
% DurationCPUTime: 0.49s
% Computational Cost: add. (262->71), mult. (477->127), div. (0->0), fcn. (514->6), ass. (0->51)
t30 = cos(pkin(9));
t24 = -t30 * pkin(5) - pkin(4);
t60 = 0.2e1 * t24;
t35 = sin(qJ(4));
t59 = -0.2e1 * t35;
t58 = 0.2e1 * t35;
t37 = cos(qJ(4));
t57 = 0.2e1 * t37;
t56 = t37 * pkin(4);
t29 = sin(pkin(9));
t34 = sin(qJ(6));
t36 = cos(qJ(6));
t16 = t34 * t29 - t36 * t30;
t55 = t16 * t35;
t17 = t36 * t29 + t34 * t30;
t54 = t17 * t35;
t28 = t37 ^ 2;
t31 = -pkin(7) + qJ(2);
t53 = t28 * t31;
t52 = t29 * t31;
t51 = t29 * t37;
t22 = t30 * t37;
t50 = t35 * t31;
t49 = t37 * t16;
t10 = t37 * t17;
t48 = t37 * t31;
t32 = pkin(1) + qJ(3);
t47 = pkin(8) + qJ(5);
t18 = t35 * pkin(4) - t37 * qJ(5) + t32;
t8 = t29 * t18 + t30 * t50;
t46 = t29 ^ 2 + t30 ^ 2;
t27 = t35 ^ 2;
t45 = -t27 - t28;
t44 = t46 * qJ(5);
t14 = t30 * t18;
t7 = -t29 * t50 + t14;
t43 = -t8 * t29 - t7 * t30;
t42 = -t7 * t29 + t8 * t30;
t41 = -qJ(5) * t35 - t56;
t40 = (qJ(2) ^ 2);
t38 = 2 * qJ(2);
t20 = t47 * t30;
t19 = t47 * t29;
t15 = (pkin(5) * t29 - t31) * t37;
t6 = -t34 * t19 + t36 * t20;
t5 = -t36 * t19 - t34 * t20;
t4 = -pkin(8) * t51 + t8;
t3 = -pkin(8) * t22 + t14 + (pkin(5) - t52) * t35;
t2 = t34 * t3 + t36 * t4;
t1 = t36 * t3 - t34 * t4;
t9 = [1, 0, 0, -2 * pkin(1), t38, pkin(1) ^ 2 + t40, t38, 0.2e1 * t32, t32 ^ 2 + t40, t28, t37 * t59, 0, 0, 0, t32 * t58, t32 * t57, -0.2e1 * t28 * t52 + 0.2e1 * t7 * t35, -0.2e1 * t30 * t53 - 0.2e1 * t8 * t35, t43 * t57, t28 * t31 ^ 2 + t7 ^ 2 + t8 ^ 2, t49 ^ 2, 0.2e1 * t49 * t10, -t49 * t58, t10 * t59, t27, 0.2e1 * t1 * t35 + 0.2e1 * t15 * t10, -0.2e1 * t15 * t49 - 0.2e1 * t2 * t35; 0, 0, 0, 1, 0, -pkin(1), 0, -1, -t32, 0, 0, 0, 0, 0, -t35, -t37, -t30 * t35, t29 * t35, t46 * t37, t43, 0, 0, 0, 0, 0, t55, t54; 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, qJ(2), 0, 0, 0, 0, 0, 0, 0, t45 * t29, t45 * t30, 0, t42 * t35 + t53, 0, 0, 0, 0, 0, -t37 * t10 - t35 * t54, t35 * t55 + t37 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t27 + t28, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t35, 0, t48, -t50, t41 * t29 + t30 * t48, -t29 * t48 + t41 * t30, t42, pkin(4) * t48 + t42 * qJ(5), -t49 * t17, -t17 * t10 + t16 * t49, t54, -t55, 0, t24 * t10 + t15 * t16 + t5 * t35, t15 * t17 - t24 * t49 - t6 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t35, t22, -t51, t46 * t35, t35 * t44 + t56, 0, 0, 0, 0, 0, -t49, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t30, -0.2e1 * pkin(4) * t29, 0.2e1 * t44, t46 * qJ(5) ^ 2 + pkin(4) ^ 2, t17 ^ 2, -0.2e1 * t17 * t16, 0, 0, 0, t16 * t60, t17 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t22, 0, -t48, 0, 0, 0, 0, 0, t10, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, 0, -pkin(4), 0, 0, 0, 0, 0, t16, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t10, t35, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
