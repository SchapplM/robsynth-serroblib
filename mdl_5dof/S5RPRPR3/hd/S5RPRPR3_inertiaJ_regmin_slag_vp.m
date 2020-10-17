% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:21
% EndTime: 2020-01-03 11:36:23
% DurationCPUTime: 0.22s
% Computational Cost: add. (174->44), mult. (357->74), div. (0->0), fcn. (314->8), ass. (0->45)
t32 = sin(pkin(9));
t54 = -0.2e1 * t32;
t53 = 0.2e1 * t32;
t35 = cos(pkin(8));
t27 = t35 * pkin(1) + pkin(2);
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t33 = sin(pkin(8));
t51 = pkin(1) * t33;
t16 = -t37 * t27 - t39 * t51;
t13 = qJ(4) - t16;
t36 = sin(qJ(5));
t15 = t39 * t27 - t37 * t51;
t34 = cos(pkin(9));
t17 = -t34 * pkin(4) - t32 * pkin(7) - pkin(3);
t4 = -t15 + t17;
t38 = cos(qJ(5));
t46 = t38 * t34;
t3 = t13 * t46 + t36 * t4;
t30 = t32 ^ 2;
t48 = t30 * t38;
t52 = t13 * t48 + t3 * t34;
t14 = -pkin(3) - t15;
t50 = pkin(3) - t14;
t28 = t30 * qJ(4);
t42 = qJ(4) * t34;
t9 = t36 * t17 + t38 * t42;
t49 = t38 * t28 + t9 * t34;
t10 = t30 * t13;
t47 = t36 * t32;
t24 = t36 * t34;
t31 = t34 ^ 2;
t45 = t31 * t13 + t10;
t44 = t31 * qJ(4) + t28;
t43 = t30 + t31;
t26 = t38 * t32;
t25 = t38 ^ 2 * t30;
t21 = t36 * t28;
t20 = -0.2e1 * t36 * t48;
t19 = t46 * t54;
t18 = t24 * t53;
t8 = t38 * t17 - t36 * t42;
t6 = t36 * t10;
t2 = -t13 * t24 + t38 * t4;
t1 = [1, 0, 0, (t33 ^ 2 + t35 ^ 2) * pkin(1) ^ 2, 1, 0.2e1 * t15, 0.2e1 * t16, -0.2e1 * t14 * t34, t14 * t53, 0.2e1 * t45, t43 * t13 ^ 2 + t14 ^ 2, t25, t20, t19, t18, t31, -0.2e1 * t2 * t34 + 0.2e1 * t6, 0.2e1 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t15, t16, t50 * t34, -t50 * t32, t44 + t45, t43 * t13 * qJ(4) - t14 * pkin(3), t25, t20, t19, t18, t31, t21 + t6 + (-t2 - t8) * t34, t49 + t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t34, pkin(3) * t54, 0.2e1 * t44, t43 * qJ(4) ^ 2 + pkin(3) ^ 2, t25, t20, t19, t18, t31, -0.2e1 * t8 * t34 + 0.2e1 * t21, 0.2e1 * t49; 0, 0, 0, 0, 0, 0, 0, -t34, t32, 0, t14, 0, 0, 0, 0, 0, -t46, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t34, t32, 0, -pkin(3), 0, 0, 0, 0, 0, -t46, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t47, -t34, t2, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t47, -t34, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t1;
