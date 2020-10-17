% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:31
% EndTime: 2019-12-31 20:01:33
% DurationCPUTime: 0.35s
% Computational Cost: add. (432->64), mult. (819->133), div. (0->0), fcn. (907->6), ass. (0->46)
t33 = cos(qJ(4));
t54 = -0.2e1 * t33;
t52 = cos(qJ(2));
t43 = t52 * pkin(6);
t21 = t52 * qJ(3) + t43;
t29 = sin(pkin(8));
t30 = cos(pkin(8));
t32 = sin(qJ(2));
t41 = (-qJ(3) - pkin(6)) * t32;
t12 = t30 * t21 + t29 * t41;
t31 = sin(qJ(4));
t18 = t29 * t32 - t30 * t52;
t19 = t29 * t52 + t30 * t32;
t26 = -t52 * pkin(2) - pkin(1);
t8 = t18 * pkin(3) - t19 * pkin(7) + t26;
t4 = t33 * t12 + t31 * t8;
t53 = t18 * pkin(4);
t24 = t29 * pkin(2) + pkin(7);
t51 = t18 * t24;
t13 = t31 * t18;
t50 = t31 * t19;
t49 = t31 * t24;
t48 = t31 * t33;
t14 = t33 * t18;
t15 = t33 * t19;
t47 = t33 * t24;
t27 = t31 ^ 2;
t28 = t33 ^ 2;
t46 = t27 + t28;
t45 = t18 * qJ(5);
t44 = 0.2e1 * t52;
t25 = -t30 * pkin(2) - pkin(3);
t42 = t31 * t12 - t33 * t8;
t10 = t29 * t21 - t30 * t41;
t1 = t45 + t4;
t2 = t42 - t53;
t40 = t1 * t33 + t2 * t31;
t39 = t1 * t31 - t2 * t33;
t38 = t33 * pkin(4) + t31 * qJ(5);
t37 = pkin(4) * t31 - t33 * qJ(5);
t16 = t25 - t38;
t36 = t16 * t19 - t51;
t35 = t19 * t25 - t51;
t17 = t19 ^ 2;
t5 = t37 * t19 + t10;
t3 = [1, 0, 0, t32 ^ 2, t32 * t44, 0, 0, 0, pkin(1) * t44, -0.2e1 * pkin(1) * t32, 0.2e1 * t10 * t19 - 0.2e1 * t12 * t18, t10 ^ 2 + t12 ^ 2 + t26 ^ 2, t28 * t17, -0.2e1 * t17 * t48, 0.2e1 * t18 * t15, -0.2e1 * t18 * t50, t18 ^ 2, 0.2e1 * t10 * t50 - 0.2e1 * t18 * t42, 0.2e1 * t10 * t15 - 0.2e1 * t4 * t18, -0.2e1 * t2 * t18 + 0.2e1 * t5 * t50, -0.2e1 * t39 * t19, 0.2e1 * t1 * t18 - 0.2e1 * t5 * t15, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t32, t52, 0, -t32 * pkin(6), -t43, (-t18 * t29 - t19 * t30) * pkin(2), (-t10 * t30 + t12 * t29) * pkin(2), t31 * t15, (-t27 + t28) * t19, t13, t14, 0, -t10 * t33 + t35 * t31, t10 * t31 + t35 * t33, t36 * t31 - t5 * t33, t40, -t5 * t31 - t36 * t33, t5 * t16 + t40 * t24; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t29 ^ 2 + t30 ^ 2) * pkin(2) ^ 2, t27, 0.2e1 * t48, 0, 0, 0, t25 * t54, 0.2e1 * t25 * t31, t16 * t54, 0.2e1 * t46 * t24, -0.2e1 * t16 * t31, t46 * t24 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, t14, -t13, t14, -t46 * t19, t13, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t50, t18, -t42, -t4, -t42 + 0.2e1 * t53, -t38 * t19, 0.2e1 * t45 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t33, 0, -t49, -t47, -t49, -t37, t47, -t37 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, t33, 0, t31, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t15, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
