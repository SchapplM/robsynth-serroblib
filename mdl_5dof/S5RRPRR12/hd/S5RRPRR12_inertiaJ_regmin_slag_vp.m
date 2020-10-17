% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:30:37
% EndTime: 2019-12-31 20:30:38
% DurationCPUTime: 0.34s
% Computational Cost: add. (212->63), mult. (409->109), div. (0->0), fcn. (461->6), ass. (0->54)
t29 = sin(qJ(4));
t30 = sin(qJ(2));
t32 = cos(qJ(4));
t33 = cos(qJ(2));
t9 = t30 * t29 + t33 * t32;
t58 = 0.2e1 * t9;
t28 = sin(qJ(5));
t57 = -0.2e1 * t28;
t56 = -0.2e1 * t30;
t31 = cos(qJ(5));
t55 = 0.2e1 * t31;
t54 = 0.2e1 * t33;
t53 = t28 * t9;
t52 = t31 * t9;
t20 = t30 * pkin(6);
t16 = -t30 * pkin(7) + t20;
t21 = t33 * pkin(6);
t17 = -t33 * pkin(7) + t21;
t5 = -t32 * t16 + t29 * t17;
t51 = t5 * t28;
t50 = t5 * t31;
t34 = -pkin(2) - pkin(3);
t13 = t29 * qJ(3) - t32 * t34;
t11 = pkin(4) + t13;
t49 = pkin(4) + t11;
t10 = -t33 * t29 + t30 * t32;
t48 = t28 * t10;
t47 = t28 * t31;
t46 = t31 * t10;
t45 = t32 * t28;
t44 = t32 * t31;
t25 = t30 ^ 2;
t43 = t33 ^ 2 + t25;
t42 = -0.2e1 * t10 * t9;
t41 = -0.2e1 * t47;
t15 = -t33 * pkin(2) - t30 * qJ(3) - pkin(1);
t40 = t28 * t46;
t7 = t33 * pkin(3) - t15;
t39 = -pkin(4) * t10 - pkin(8) * t9;
t14 = t32 * qJ(3) + t29 * t34;
t12 = -pkin(8) + t14;
t38 = t10 * t11 - t12 * t9;
t37 = -t10 * t32 - t29 * t9;
t36 = -t30 * pkin(2) + t33 * qJ(3);
t26 = t31 ^ 2;
t24 = t28 ^ 2;
t18 = 0.2e1 * t47;
t8 = t10 ^ 2;
t6 = t29 * t16 + t32 * t17;
t4 = (t24 - t26) * t10;
t3 = t9 * pkin(4) - t10 * pkin(8) + t7;
t2 = t28 * t3 + t31 * t6;
t1 = -t28 * t6 + t31 * t3;
t19 = [1, 0, 0, t25, t30 * t54, 0, 0, 0, pkin(1) * t54, pkin(1) * t56, -0.2e1 * t15 * t33, 0.2e1 * t43 * pkin(6), t15 * t56, t43 * pkin(6) ^ 2 + t15 ^ 2, t8, t42, 0, 0, 0, t7 * t58, 0.2e1 * t7 * t10, t26 * t8, t8 * t41, t46 * t58, t28 * t42, t9 ^ 2, 0.2e1 * t1 * t9 + 0.2e1 * t5 * t48, -0.2e1 * t2 * t9 + 0.2e1 * t5 * t46; 0, 0, 0, 0, 0, t30, t33, 0, -t20, -t21, -t20, t36, t21, t36 * pkin(6), 0, 0, -t10, t9, 0, t5, t6, -t40, t4, -t53, -t52, 0, t38 * t28 + t50, t38 * t31 - t51; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t13, 0.2e1 * t14, t24, t18, 0, 0, 0, t11 * t55, t11 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t28, t37 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, -t32, t29, 0, 0, 0, 0, 0, -t44, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, -t5, -t6, t40, -t4, t53, t52, 0, t39 * t28 - t50, t39 * t31 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t13, -t14, -t24, t41, 0, 0, 0, -t49 * t31, t49 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t29, 0, 0, 0, 0, 0, t44, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t24, t18, 0, 0, 0, pkin(4) * t55, pkin(4) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t48, t9, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t31, 0, -t28 * t12, -t31 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t29, -t31 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t31, 0, -t28 * pkin(8), -t31 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t19;
