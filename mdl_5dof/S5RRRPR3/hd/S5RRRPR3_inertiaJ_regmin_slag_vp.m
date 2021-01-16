% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:59:36
% EndTime: 2021-01-15 22:59:38
% DurationCPUTime: 0.33s
% Computational Cost: add. (401->57), mult. (760->107), div. (0->0), fcn. (874->8), ass. (0->60)
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t48 = sin(qJ(3));
t51 = cos(qJ(3));
t31 = t45 * t48 - t46 * t51;
t42 = -t51 * pkin(3) - pkin(2);
t23 = t31 * pkin(4) + t42;
t59 = cos(qJ(2)) * pkin(1);
t22 = t23 - t59;
t70 = 0.2e1 * t22;
t69 = 0.2e1 * t23;
t34 = t42 - t59;
t68 = 0.2e1 * t34;
t67 = 0.2e1 * t42;
t66 = 0.2e1 * t51;
t61 = sin(qJ(2)) * pkin(1);
t40 = pkin(7) + t61;
t29 = (-qJ(4) - t40) * t48;
t43 = t51 * qJ(4);
t57 = t51 * t40;
t30 = t43 + t57;
t14 = t46 * t29 - t45 * t30;
t15 = t45 * t29 + t46 * t30;
t32 = t45 * t51 + t46 * t48;
t65 = -t14 * t32 - t15 * t31;
t64 = t32 * pkin(8);
t63 = t45 * pkin(3);
t62 = t46 * pkin(3);
t60 = t51 * pkin(7);
t41 = -pkin(2) - t59;
t58 = pkin(2) - t41;
t35 = (-qJ(4) - pkin(7)) * t48;
t36 = t43 + t60;
t20 = t46 * t35 - t45 * t36;
t21 = t45 * t35 + t46 * t36;
t56 = -t20 * t32 - t21 * t31;
t55 = t22 + t23;
t54 = t34 + t42;
t50 = cos(qJ(5));
t47 = sin(qJ(5));
t44 = t48 ^ 2;
t39 = pkin(4) + t62;
t37 = t48 * t66;
t28 = t31 * pkin(8);
t26 = -t47 * t39 - t50 * t63;
t25 = t50 * t39 - t47 * t63;
t19 = -t47 * t31 + t50 * t32;
t18 = t50 * t31 + t47 * t32;
t17 = t19 ^ 2;
t16 = (-t31 * t45 - t32 * t46) * pkin(3);
t11 = t21 - t28;
t10 = t20 - t64;
t7 = t15 - t28;
t6 = t14 - t64;
t5 = -0.2e1 * t19 * t18;
t4 = -t47 * t10 - t50 * t11;
t3 = t50 * t10 - t47 * t11;
t2 = -t47 * t6 - t50 * t7;
t1 = -t47 * t7 + t50 * t6;
t8 = [1, 0, 0, 1, 0.2e1 * t59, -0.2e1 * t61, t44, t37, 0, 0, 0, -0.2e1 * t41 * t51, 0.2e1 * t41 * t48, t31 * t68, t32 * t68, 0.2e1 * t65, t14 ^ 2 + t15 ^ 2 + t34 ^ 2, t17, t5, 0, 0, 0, t18 * t70, t19 * t70; 0, 0, 0, 1, t59, -t61, t44, t37, 0, 0, 0, t58 * t51, -t58 * t48, t54 * t31, t54 * t32, t56 + t65, t14 * t20 + t15 * t21 + t34 * t42, t17, t5, 0, 0, 0, t55 * t18, t55 * t19; 0, 0, 0, 1, 0, 0, t44, t37, 0, 0, 0, pkin(2) * t66, -0.2e1 * pkin(2) * t48, t31 * t67, t32 * t67, 0.2e1 * t56, t20 ^ 2 + t21 ^ 2 + t42 ^ 2, t17, t5, 0, 0, 0, t18 * t69, t19 * t69; 0, 0, 0, 0, 0, 0, 0, 0, t48, t51, 0, -t48 * t40, -t57, t14, -t15, t16, (t14 * t46 + t15 * t45) * pkin(3), 0, 0, t19, -t18, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, t48, t51, 0, -t48 * pkin(7), -t60, t20, -t21, t16, (t20 * t46 + t21 * t45) * pkin(3), 0, 0, t19, -t18, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t62, -0.2e1 * t63, 0, (t45 ^ 2 + t46 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t25, 0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t32, 0, t34, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t32, 0, t42, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t25, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
