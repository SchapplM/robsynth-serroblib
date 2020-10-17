% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:34:22
% EndTime: 2019-05-05 17:34:24
% DurationCPUTime: 0.50s
% Computational Cost: add. (647->83), mult. (1139->149), div. (0->0), fcn. (1284->8), ass. (0->57)
t36 = sin(pkin(10));
t38 = cos(pkin(10));
t41 = sin(qJ(3));
t63 = cos(qJ(3));
t26 = t36 * t63 + t38 * t41;
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t47 = pkin(5) * t40 - t42 * qJ(6);
t67 = t47 * t26;
t24 = t36 * t41 - t38 * t63;
t22 = t24 ^ 2;
t66 = 0.2e1 * t41;
t65 = -0.2e1 * t42;
t37 = sin(pkin(9));
t31 = t37 * pkin(1) + pkin(7);
t54 = t63 * t31;
t21 = t63 * qJ(4) + t54;
t51 = (-qJ(4) - t31) * t41;
t12 = t38 * t21 + t36 * t51;
t39 = cos(pkin(9));
t33 = -t39 * pkin(1) - pkin(2);
t27 = -t63 * pkin(3) + t33;
t9 = t24 * pkin(4) - t26 * pkin(8) + t27;
t4 = t42 * t12 + t40 * t9;
t64 = t24 * pkin(5);
t30 = t36 * pkin(3) + pkin(8);
t62 = t24 * t30;
t34 = t40 ^ 2;
t61 = t34 * t26;
t15 = t40 * t24;
t60 = t40 * t26;
t59 = t40 * t30;
t58 = t40 * t42;
t18 = t42 * t24;
t19 = t42 * t26;
t57 = t42 * t30;
t35 = t42 ^ 2;
t56 = t34 + t35;
t55 = t24 * qJ(6);
t32 = -t38 * pkin(3) - pkin(4);
t53 = t40 * t12 - t42 * t9;
t52 = t56 * t30;
t10 = t36 * t21 - t38 * t51;
t1 = t55 + t4;
t2 = t53 - t64;
t50 = t1 * t42 + t2 * t40;
t49 = t1 * t40 - t2 * t42;
t48 = t42 * pkin(5) + t40 * qJ(6);
t20 = t32 - t48;
t46 = t20 * t26 - t62;
t45 = t26 * t32 - t62;
t23 = t26 ^ 2;
t17 = t35 * t26;
t16 = t35 * t23;
t13 = -t17 - t61;
t5 = t10 + t67;
t3 = [1, 0, 0 (t37 ^ 2 + t39 ^ 2) * pkin(1) ^ 2, t41 ^ 2, t63 * t66, 0, 0, 0, -0.2e1 * t33 * t63, t33 * t66, 0.2e1 * t10 * t26 - 0.2e1 * t12 * t24, t10 ^ 2 + t12 ^ 2 + t27 ^ 2, t16, -0.2e1 * t23 * t58, 0.2e1 * t24 * t19, -0.2e1 * t24 * t60, t22, 0.2e1 * t10 * t60 - 0.2e1 * t24 * t53, 0.2e1 * t10 * t19 - 0.2e1 * t4 * t24, -0.2e1 * t2 * t24 + 0.2e1 * t5 * t60, -0.2e1 * t49 * t26, 0.2e1 * t1 * t24 - 0.2e1 * t5 * t19, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t24 + t12 * t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t24 + t50 * t26; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t23 + t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t23 + t16 + t22; 0, 0, 0, 0, 0, 0, t41, t63, 0, -t41 * t31, -t54 (-t24 * t36 - t26 * t38) * pkin(3) (-t10 * t38 + t12 * t36) * pkin(3), t40 * t19, t17 - t61, t15, t18, 0, -t10 * t42 + t45 * t40, t10 * t40 + t45 * t42, t46 * t40 - t5 * t42, t50, -t5 * t40 - t46 * t42, t5 * t20 + t50 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t41, 0 (-t24 * t38 + t26 * t36) * pkin(3), 0, 0, 0, 0, 0, -t18, t15, -t18, -t13, -t15, t24 * t20 + t26 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t36 ^ 2 + t38 ^ 2) * pkin(3) ^ 2, t34, 0.2e1 * t58, 0, 0, 0, t32 * t65, 0.2e1 * t32 * t40, t20 * t65, 0.2e1 * t52, -0.2e1 * t20 * t40, t56 * t30 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, t18, -t15, t18, t13, t15, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t60, t24, -t53, -t4, -t53 + 0.2e1 * t64, -t48 * t26, 0.2e1 * t55 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t19, -t60, 0, t19, -t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t42, 0, -t59, -t57, -t59, -t47, t57, -t47 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t40, t42, 0, t40, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t19, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
