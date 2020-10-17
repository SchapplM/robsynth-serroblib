% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:38:29
% EndTime: 2019-05-05 14:38:29
% DurationCPUTime: 0.38s
% Computational Cost: add. (290->56), mult. (498->93), div. (0->0), fcn. (592->6), ass. (0->45)
t33 = sin(pkin(9));
t34 = cos(pkin(9));
t37 = sin(qJ(4));
t51 = cos(qJ(4));
t19 = t51 * t33 + t37 * t34;
t15 = t19 ^ 2;
t18 = t37 * t33 - t51 * t34;
t56 = t18 ^ 2;
t57 = t56 + t15;
t55 = 2 * qJ(2);
t54 = 2 * qJ(5);
t53 = pkin(4) + pkin(8);
t35 = -pkin(1) - qJ(3);
t52 = -pkin(7) + t35;
t50 = t18 * t19;
t36 = sin(qJ(6));
t49 = t36 * t18;
t10 = t36 * t19;
t38 = cos(qJ(6));
t48 = t38 * t18;
t12 = t38 * t19;
t47 = t38 * t36;
t24 = t33 ^ 2 + t34 ^ 2;
t46 = qJ(5) * t19;
t25 = t33 * pkin(3) + qJ(2);
t45 = -0.2e1 * t50;
t22 = t52 * t33;
t23 = t52 * t34;
t7 = t37 * t22 - t51 * t23;
t8 = t51 * t22 + t37 * t23;
t44 = t7 * t18 + t8 * t19;
t43 = t18 * qJ(5) + t25;
t42 = t18 * pkin(4) - t46;
t41 = -t18 * t53 + t46;
t40 = qJ(2) ^ 2;
t32 = t38 ^ 2;
t31 = t36 ^ 2;
t17 = t24 * t35;
t6 = t19 * pkin(4) + t43;
t5 = -t19 * pkin(5) + t8;
t4 = -t18 * pkin(5) + t7;
t3 = t53 * t19 + t43;
t2 = t38 * t3 + t36 * t4;
t1 = -t36 * t3 + t38 * t4;
t9 = [1, 0, 0, -2 * pkin(1), t55, pkin(1) ^ 2 + t40, t33 * t55, t34 * t55, -0.2e1 * t17, t24 * t35 ^ 2 + t40, t56, 0.2e1 * t50, 0, 0, 0, 0.2e1 * t25 * t19, -0.2e1 * t25 * t18, -0.2e1 * t44, -0.2e1 * t6 * t19, 0.2e1 * t6 * t18, t6 ^ 2 + t7 ^ 2 + t8 ^ 2, t31 * t15, 0.2e1 * t15 * t47, t36 * t45, t38 * t45, t56, -0.2e1 * t1 * t18 - 0.2e1 * t5 * t12, 0.2e1 * t5 * t10 + 0.2e1 * t2 * t18; 0, 0, 0, 1, 0, -pkin(1), 0, 0, -t24, t17, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, t44, 0, 0, 0, 0, 0, -t57 * t38, t57 * t36; 0, 0, 0, 0, 0, 1, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t33, t34, 0, qJ(2), 0, 0, 0, 0, 0, t19, -t18, 0, -t19, t18, t6, 0, 0, 0, 0, 0, t49, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, -t7, -t8, t42, t7, t8, -t7 * pkin(4) + t8 * qJ(5), t19 * t47 (-t31 + t32) * t19, -t48, t49, 0, t5 * t36 - t41 * t38, t41 * t36 + t5 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, t18, t19, -t42, 0, 0, 0, 0, 0, t10, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t54, pkin(4) ^ 2 + (qJ(5) ^ 2) t32, -0.2e1 * t47, 0, 0, 0, t36 * t54, t38 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, t7, 0, 0, 0, 0, 0, -t48, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t12, -t18, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t36, 0, -t38 * t53, t36 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
