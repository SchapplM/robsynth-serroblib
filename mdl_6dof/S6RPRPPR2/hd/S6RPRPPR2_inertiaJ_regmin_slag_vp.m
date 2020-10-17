% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:34:12
% EndTime: 2019-05-05 16:34:13
% DurationCPUTime: 0.38s
% Computational Cost: add. (353->56), mult. (633->107), div. (0->0), fcn. (728->8), ass. (0->48)
t31 = sin(pkin(10));
t33 = cos(pkin(10));
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t17 = t31 * t36 - t33 * t38;
t15 = t17 ^ 2;
t19 = t31 * t38 + t33 * t36;
t16 = t19 ^ 2;
t56 = t15 + t16;
t34 = cos(pkin(9));
t28 = -t34 * pkin(1) - pkin(2);
t21 = -t38 * pkin(3) + t28;
t41 = -t19 * qJ(5) + t21;
t6 = t17 * pkin(4) + t41;
t55 = -0.2e1 * t6;
t23 = t31 * pkin(3) + qJ(5);
t54 = 0.2e1 * t23;
t53 = 0.2e1 * t36;
t52 = t23 * t17;
t35 = sin(qJ(6));
t51 = t35 * t17;
t50 = t35 * t19;
t37 = cos(qJ(6));
t12 = t37 * t17;
t13 = t37 * t19;
t49 = t37 * t35;
t32 = sin(pkin(9));
t26 = t32 * pkin(1) + pkin(7);
t48 = qJ(4) + t26;
t47 = 0.2e1 * t17 * t19;
t14 = t48 * t38;
t45 = t48 * t36;
t7 = t31 * t14 + t33 * t45;
t9 = t33 * t14 - t31 * t45;
t46 = t7 ^ 2 + t9 ^ 2;
t27 = -t33 * pkin(3) - pkin(4);
t44 = t7 * t17 + t9 * t19;
t22 = -pkin(8) + t27;
t43 = -t19 * t22 + t52;
t42 = -0.2e1 * t9 * t17 + 0.2e1 * t7 * t19;
t30 = t37 ^ 2;
t29 = t35 ^ 2;
t5 = -t17 * pkin(5) + t9;
t4 = t19 * pkin(5) + t7;
t3 = (pkin(4) + pkin(8)) * t17 + t41;
t2 = t37 * t3 + t35 * t4;
t1 = -t35 * t3 + t37 * t4;
t8 = [1, 0, 0 (t32 ^ 2 + t34 ^ 2) * pkin(1) ^ 2, t36 ^ 2, t38 * t53, 0, 0, 0, -0.2e1 * t28 * t38, t28 * t53, t42, t21 ^ 2 + t46, t42, t17 * t55, t19 * t55, t6 ^ 2 + t46, t29 * t15, 0.2e1 * t15 * t49, t35 * t47, t37 * t47, t16, 0.2e1 * t1 * t19 - 0.2e1 * t5 * t12, -0.2e1 * t2 * t19 + 0.2e1 * t5 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, 0, 0, t56, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t36, t38, 0, -t36 * t26, -t38 * t26 (-t17 * t31 - t19 * t33) * pkin(3) (t31 * t9 - t33 * t7) * pkin(3), t27 * t19 - t52, t7, t9, t9 * t23 + t7 * t27, t17 * t49 (-t29 + t30) * t17, t13, -t50, 0, t5 * t35 - t37 * t43, t35 * t43 + t5 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t36, 0 (-t17 * t33 + t19 * t31) * pkin(3), 0, t17, t19, t17 * t27 + t19 * t23, 0, 0, 0, 0, 0, t50, t13; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t31 ^ 2 + t33 ^ 2) * pkin(3) ^ 2, 0, 0.2e1 * t27, t54, t23 ^ 2 + t27 ^ 2, t30, -0.2e1 * t49, 0, 0, 0, t35 * t54, t37 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t17, -t19, t6, 0, 0, 0, 0, 0, -t50, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, t7, 0, 0, 0, 0, 0, t13, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t27, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t12, t19, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t35, 0, t37 * t22, -t35 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
