% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR15_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:24
% EndTime: 2019-12-31 18:37:25
% DurationCPUTime: 0.33s
% Computational Cost: add. (215->64), mult. (442->121), div. (0->0), fcn. (472->6), ass. (0->46)
t32 = sin(qJ(3));
t53 = -0.2e1 * t32;
t30 = cos(pkin(8));
t24 = -pkin(4) * t30 - pkin(3);
t54 = 0.2e1 * t24;
t52 = 2 * qJ(2);
t34 = cos(qJ(3));
t51 = pkin(3) * t34;
t29 = sin(pkin(8));
t31 = sin(qJ(5));
t33 = cos(qJ(5));
t16 = t29 * t31 - t33 * t30;
t50 = t16 * t32;
t17 = t29 * t33 + t30 * t31;
t49 = t17 * t32;
t28 = t34 ^ 2;
t35 = -pkin(1) - pkin(6);
t48 = t28 * t35;
t47 = t29 * t34;
t46 = t29 * t35;
t22 = t30 * t34;
t45 = t32 * t35;
t44 = t34 * t16;
t10 = t34 * t17;
t43 = t34 * t35;
t42 = pkin(7) + qJ(4);
t18 = pkin(3) * t32 - qJ(4) * t34 + qJ(2);
t8 = t29 * t18 + t30 * t45;
t41 = t29 ^ 2 + t30 ^ 2;
t27 = t32 ^ 2;
t40 = -t27 - t28;
t39 = t41 * qJ(4);
t14 = t30 * t18;
t7 = -t29 * t45 + t14;
t38 = -t29 * t7 + t30 * t8;
t37 = -qJ(4) * t32 - t51;
t20 = t42 * t30;
t19 = t42 * t29;
t15 = (pkin(4) * t29 - t35) * t34;
t6 = -t19 * t31 + t20 * t33;
t5 = -t19 * t33 - t20 * t31;
t4 = -pkin(7) * t47 + t8;
t3 = -pkin(7) * t22 + t14 + (pkin(4) - t46) * t32;
t2 = t3 * t31 + t33 * t4;
t1 = t3 * t33 - t31 * t4;
t9 = [1, 0, 0, -2 * pkin(1), t52, pkin(1) ^ 2 + qJ(2) ^ 2, t28, t34 * t53, 0, 0, 0, t32 * t52, t34 * t52, -0.2e1 * t28 * t46 + 0.2e1 * t32 * t7, -0.2e1 * t30 * t48 - 0.2e1 * t32 * t8, 0.2e1 * (-t29 * t8 - t30 * t7) * t34, t28 * t35 ^ 2 + t7 ^ 2 + t8 ^ 2, t44 ^ 2, 0.2e1 * t44 * t10, t44 * t53, t10 * t53, t27, 0.2e1 * t1 * t32 + 0.2e1 * t10 * t15, -0.2e1 * t15 * t44 - 0.2e1 * t2 * t32; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, t40 * t29, t40 * t30, 0, t38 * t32 + t48, 0, 0, 0, 0, 0, -t10 * t34 - t32 * t49, t32 * t50 + t34 * t44; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t27 + t28, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, t43, -t45, t37 * t29 + t30 * t43, -t29 * t43 + t37 * t30, t38, pkin(3) * t43 + t38 * qJ(4), -t44 * t17, -t10 * t17 + t16 * t44, t49, -t50, 0, t10 * t24 + t15 * t16 + t32 * t5, t15 * t17 - t24 * t44 - t32 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, t22, -t47, t41 * t32, t32 * t39 + t51, 0, 0, 0, 0, 0, -t44, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t30, -0.2e1 * pkin(3) * t29, 0.2e1 * t39, t41 * qJ(4) ^ 2 + pkin(3) ^ 2, t17 ^ 2, -0.2e1 * t17 * t16, 0, 0, 0, t16 * t54, t17 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t22, 0, -t43, 0, 0, 0, 0, 0, t10, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t29, 0, -pkin(3), 0, 0, 0, 0, 0, t16, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t10, t32, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
