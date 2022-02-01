% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:02:06
% EndTime: 2022-01-20 12:02:08
% DurationCPUTime: 0.36s
% Computational Cost: add. (225->55), mult. (434->75), div. (0->0), fcn. (478->8), ass. (0->58)
t45 = sin(qJ(4));
t63 = pkin(3) * t45;
t44 = sin(qJ(5));
t62 = t44 * pkin(4);
t46 = sin(qJ(3));
t61 = t46 * pkin(2);
t60 = sin(qJ(2)) * pkin(1);
t48 = cos(qJ(5));
t59 = t48 * pkin(4);
t49 = cos(qJ(4));
t58 = t49 * pkin(8);
t41 = cos(qJ(2)) * pkin(1);
t37 = t41 + pkin(2);
t50 = cos(qJ(3));
t53 = -t50 * t37 + t46 * t60;
t18 = -pkin(3) + t53;
t57 = t18 * t49;
t40 = t50 * pkin(2);
t36 = -t40 - pkin(3);
t56 = t36 * t49;
t52 = t50 * t60;
t21 = -t46 * t37 - t52;
t19 = pkin(8) - t21;
t55 = t49 * t19;
t35 = pkin(8) + t61;
t54 = t49 * t35;
t38 = -t49 * pkin(4) - pkin(3);
t43 = t45 ^ 2;
t42 = pkin(3) * t49;
t39 = t49 * pkin(9);
t33 = 0.2e1 * t45 * t49;
t31 = t36 * t45;
t29 = t39 + t58;
t28 = (-pkin(8) - pkin(9)) * t45;
t27 = t38 - t40;
t26 = t44 * t49 + t48 * t45;
t25 = t44 * t45 - t48 * t49;
t24 = t26 ^ 2;
t23 = t39 + t54;
t22 = (-pkin(9) - t35) * t45;
t17 = t38 * t26;
t16 = t38 * t25;
t15 = t18 * t45;
t14 = t38 + t53;
t13 = t27 * t26;
t12 = t27 * t25;
t11 = -t44 * t28 - t48 * t29;
t10 = t48 * t28 - t44 * t29;
t9 = t39 + t55;
t8 = (-pkin(9) - t19) * t45;
t7 = -0.2e1 * t26 * t25;
t6 = t14 * t26;
t5 = t14 * t25;
t4 = -t44 * t22 - t48 * t23;
t3 = t48 * t22 - t44 * t23;
t2 = -t44 * t8 - t48 * t9;
t1 = -t44 * t9 + t48 * t8;
t20 = [1, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t60, 1, -0.2e1 * t53, 0.2e1 * t21, t43, t33, 0, 0, 0, -0.2e1 * t57, 0.2e1 * t15, t24, t7, 0, 0, 0, 0.2e1 * t5, 0.2e1 * t6; 0, 0, 0, 1, t41, -t60, 1, t40 - t53, -t52 + (-pkin(2) - t37) * t46, t43, t33, 0, 0, 0, (-t18 - t36) * t49, t31 + t15, t24, t7, 0, 0, 0, t12 + t5, t13 + t6; 0, 0, 0, 1, 0, 0, 1, 0.2e1 * t40, -0.2e1 * t61, t43, t33, 0, 0, 0, -0.2e1 * t56, 0.2e1 * t31, t24, t7, 0, 0, 0, 0.2e1 * t12, 0.2e1 * t13; 0, 0, 0, 0, 0, 0, 1, -t53, t21, t43, t33, 0, 0, 0, t42 - t57, t15 - t63, t24, t7, 0, 0, 0, t16 + t5, t17 + t6; 0, 0, 0, 0, 0, 0, 1, t40, -t61, t43, t33, 0, 0, 0, t42 - t56, t31 - t63, t24, t7, 0, 0, 0, t16 + t12, t17 + t13; 0, 0, 0, 0, 0, 0, 1, 0, 0, t43, t33, 0, 0, 0, 0.2e1 * t42, -0.2e1 * t63, t24, t7, 0, 0, 0, 0.2e1 * t16, 0.2e1 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t49, 0, -t45 * t19, -t55, 0, 0, t26, -t25, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t49, 0, -t45 * t35, -t54, 0, 0, t26, -t25, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t49, 0, -t45 * pkin(8), -t58, 0, 0, t26, -t25, 0, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t59, -0.2e1 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t10, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t59, -t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t20;
