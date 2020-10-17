% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:07:24
% EndTime: 2019-05-05 15:07:26
% DurationCPUTime: 0.56s
% Computational Cost: add. (563->77), mult. (964->129), div. (0->0), fcn. (1114->6), ass. (0->57)
t34 = sin(pkin(9));
t35 = cos(pkin(9));
t38 = sin(qJ(4));
t56 = cos(qJ(4));
t20 = t56 * t34 + t38 * t35;
t16 = t20 ^ 2;
t19 = t38 * t34 - t56 * t35;
t17 = t19 ^ 2;
t69 = t16 + t17;
t39 = cos(qJ(5));
t68 = t69 * t39;
t37 = sin(qJ(5));
t66 = t69 * t37;
t65 = -0.2e1 * t19;
t64 = 0.2e1 * t20;
t63 = -0.2e1 * t37;
t62 = 2 * qJ(2);
t36 = -pkin(1) - qJ(3);
t57 = -pkin(7) + t36;
t22 = t57 * t34;
t23 = t57 * t35;
t11 = t56 * t22 + t38 * t23;
t26 = t34 * pkin(3) + qJ(2);
t8 = t20 * pkin(4) + t19 * pkin(8) + t26;
t4 = t39 * t11 + t37 * t8;
t61 = pkin(8) * t20;
t60 = t20 * pkin(5);
t59 = t37 * pkin(8);
t58 = t39 * pkin(8);
t43 = t39 * pkin(5) + t37 * qJ(6);
t24 = -pkin(4) - t43;
t55 = t19 * t24;
t54 = t19 * t37;
t15 = t19 * t39;
t13 = t37 * t20;
t53 = t37 * t39;
t14 = t39 * t20;
t25 = t34 ^ 2 + t35 ^ 2;
t32 = t37 ^ 2;
t33 = t39 ^ 2;
t52 = t32 + t33;
t51 = t20 * qJ(6);
t50 = t19 * t64;
t49 = t37 * t11 - t39 * t8;
t48 = t52 * t20;
t47 = pkin(4) * t19 - t61;
t1 = t51 + t4;
t2 = t49 - t60;
t46 = t1 * t39 + t2 * t37;
t45 = t1 * t37 - t2 * t39;
t44 = t55 + t61;
t42 = -pkin(5) * t37 + t39 * qJ(6);
t10 = t38 * t22 - t56 * t23;
t40 = qJ(2) ^ 2;
t18 = t25 * t36;
t5 = t42 * t19 + t10;
t3 = [1, 0, 0, -2 * pkin(1), t62, pkin(1) ^ 2 + t40, t34 * t62, t35 * t62, -0.2e1 * t18, t25 * t36 ^ 2 + t40, t17, t50, 0, 0, 0, t26 * t64, t26 * t65, t33 * t17, -0.2e1 * t17 * t53, t14 * t65, t37 * t50, t16, -0.2e1 * t10 * t54 - 0.2e1 * t20 * t49, -0.2e1 * t10 * t15 - 0.2e1 * t4 * t20, -0.2e1 * t2 * t20 - 0.2e1 * t5 * t54, 0.2e1 * t45 * t19, 0.2e1 * t1 * t20 + 0.2e1 * t5 * t15, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, -t25, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t68, -t66, 0, t68, t5 * t19 + t46 * t20; 0, 0, 0, 0, 0, 1, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t16 + t17; 0, 0, 0, 0, 0, 0, t34, t35, 0, qJ(2), 0, 0, 0, 0, 0, t20, -t19, 0, 0, 0, 0, 0, t14, -t13, t14, t52 * t19, t13, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, -t10, -t11, -t19 * t53 (t32 - t33) * t19, t13, t14, 0, -t10 * t39 + t47 * t37, t10 * t37 + t47 * t39, -t44 * t37 - t5 * t39, t46, -t5 * t37 + t44 * t39, t46 * pkin(8) + t5 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, 0, 0, 0, 0, -t15, t54, -t15, t48, -t54, pkin(8) * t48 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t32, 0.2e1 * t53, 0, 0, 0, 0.2e1 * pkin(4) * t39, pkin(4) * t63, -0.2e1 * t24 * t39, 0.2e1 * t52 * pkin(8), t24 * t63, t52 * pkin(8) ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t54, t20, -t49, -t4, -t49 + 0.2e1 * t60, t43 * t19, 0.2e1 * t51 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, -t13, 0, t14, t42 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, t39, 0, t37, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t39, 0, -t59, -t58, -t59, t42, t58, t42 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t15, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
