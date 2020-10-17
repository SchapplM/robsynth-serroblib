% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:15:01
% EndTime: 2019-05-05 18:15:02
% DurationCPUTime: 0.45s
% Computational Cost: add. (586->70), mult. (1096->128), div. (0->0), fcn. (1353->10), ass. (0->56)
t38 = sin(pkin(11));
t40 = cos(pkin(11));
t44 = sin(qJ(3));
t46 = cos(qJ(3));
t29 = -t38 * t44 + t40 * t46;
t41 = cos(pkin(10));
t35 = -t41 * pkin(1) - pkin(2);
t31 = -t46 * pkin(3) + t35;
t18 = -t29 * pkin(4) + t31;
t62 = 0.2e1 * t18;
t61 = 0.2e1 * t44;
t60 = pkin(3) * t38;
t43 = sin(qJ(5));
t39 = sin(pkin(10));
t33 = t39 * pkin(1) + pkin(7);
t53 = qJ(4) + t33;
t27 = t53 * t44;
t28 = t53 * t46;
t11 = -t40 * t27 - t38 * t28;
t30 = t38 * t46 + t40 * t44;
t49 = -t30 * pkin(8) + t11;
t57 = cos(qJ(5));
t12 = -t38 * t27 + t40 * t28;
t9 = t29 * pkin(8) + t12;
t4 = t43 * t9 - t57 * t49;
t45 = cos(qJ(6));
t59 = t4 * t45;
t34 = t40 * pkin(3) + pkin(4);
t24 = t57 * t34 - t43 * t60;
t22 = -pkin(5) - t24;
t58 = pkin(5) - t22;
t16 = -t57 * t29 + t43 * t30;
t42 = sin(qJ(6));
t13 = t42 * t16;
t17 = t43 * t29 + t57 * t30;
t56 = t42 * t17;
t55 = t42 * t45;
t54 = t45 * t17;
t52 = -0.2e1 * t17 * t16;
t51 = -pkin(5) * t17 - pkin(9) * t16;
t25 = -t43 * t34 - t57 * t60;
t23 = pkin(9) - t25;
t50 = -t16 * t23 + t17 * t22;
t37 = t45 ^ 2;
t36 = t42 ^ 2;
t32 = 0.2e1 * t55;
t15 = t17 ^ 2;
t14 = t45 * t16;
t10 = t42 * t54;
t7 = (-t36 + t37) * t17;
t6 = t16 * pkin(5) - t17 * pkin(9) + t18;
t5 = t43 * t49 + t57 * t9;
t3 = t4 * t42;
t2 = t42 * t6 + t45 * t5;
t1 = -t42 * t5 + t45 * t6;
t8 = [1, 0, 0 (t39 ^ 2 + t41 ^ 2) * pkin(1) ^ 2, t44 ^ 2, t46 * t61, 0, 0, 0, -0.2e1 * t35 * t46, t35 * t61, -0.2e1 * t11 * t30 + 0.2e1 * t12 * t29, t11 ^ 2 + t12 ^ 2 + t31 ^ 2, t15, t52, 0, 0, 0, t16 * t62, t17 * t62, t37 * t15, -0.2e1 * t15 * t55, 0.2e1 * t16 * t54, t42 * t52, t16 ^ 2, 0.2e1 * t1 * t16 + 0.2e1 * t4 * t56, -0.2e1 * t2 * t16 + 0.2e1 * t4 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t29 + t12 * t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t29 ^ 2 + t30 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t44, t46, 0, -t44 * t33, -t46 * t33 (t29 * t38 - t30 * t40) * pkin(3) (t11 * t40 + t12 * t38) * pkin(3), 0, 0, t17, -t16, 0, -t4, -t5, t10, t7, t13, t14, 0, t50 * t42 - t59, t50 * t45 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t44, 0 (t29 * t40 + t30 * t38) * pkin(3), 0, 0, 0, 0, 0, -t16, -t17, 0, 0, 0, 0, 0, -t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t38 ^ 2 + t40 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t24, 0.2e1 * t25, t36, t32, 0, 0, 0, -0.2e1 * t22 * t45, 0.2e1 * t22 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, t14, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, -t4, -t5, t10, t7, t13, t14, 0, t51 * t42 - t59, t51 * t45 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, 0, 0, 0, 0, -t14, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t24, t25, t36, t32, 0, 0, 0, t58 * t45, -t58 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t36, t32, 0, 0, 0, 0.2e1 * pkin(5) * t45, -0.2e1 * pkin(5) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t56, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t45, 0, -t42 * t23, -t45 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t45, 0, -t42 * pkin(9), -t45 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
