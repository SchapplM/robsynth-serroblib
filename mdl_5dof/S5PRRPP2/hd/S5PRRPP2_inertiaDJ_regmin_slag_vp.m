% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:09
% EndTime: 2019-12-05 16:10:11
% DurationCPUTime: 0.31s
% Computational Cost: add. (347->64), mult. (949->126), div. (0->0), fcn. (825->6), ass. (0->47)
t39 = sin(qJ(2));
t40 = cos(qJ(3));
t41 = cos(qJ(2));
t57 = t41 * qJD(2);
t38 = sin(qJ(3));
t59 = t38 * qJD(3);
t64 = t39 * t59 - t40 * t57;
t37 = sin(pkin(8));
t60 = cos(pkin(8));
t51 = t60 * t38;
t24 = t37 * t40 + t51;
t21 = t24 * qJD(3);
t63 = 2 * qJD(5);
t62 = t37 * t38;
t61 = -qJ(4) - pkin(6);
t36 = t39 * qJD(2);
t58 = t40 * qJD(3);
t56 = -0.2e1 * pkin(2) * qJD(3);
t35 = pkin(3) * t59;
t54 = t41 * t59;
t34 = -t40 * pkin(3) - pkin(2);
t49 = qJD(3) * t61;
t20 = t40 * qJD(4) + t38 * t49;
t44 = -t38 * qJD(4) + t40 * t49;
t11 = t37 * t20 - t60 * t44;
t12 = t60 * t20 + t37 * t44;
t29 = t61 * t40;
t15 = -t37 * t29 - t61 * t51;
t16 = -t60 * t29 + t61 * t62;
t52 = t15 * t11 + t16 * t12;
t50 = t60 * t40;
t48 = qJD(3) * t50;
t45 = t50 - t62;
t10 = -t39 * t21 + t45 * t57;
t18 = t24 * t39;
t19 = t45 * t39;
t9 = t64 * t37 - t39 * t48 - t51 * t57;
t47 = t10 * t16 + t18 * t11 + t19 * t12 - t9 * t15;
t22 = -t37 * t59 + t48;
t46 = t10 * t45 + t18 * t22 - t19 * t21 - t9 * t24;
t43 = 0.2e1 * t11 * t24 + 0.2e1 * t12 * t45 + 0.2e1 * t15 * t22 - 0.2e1 * t16 * t21;
t42 = 0.2e1 * t19 * t10 - 0.2e1 * t18 * t9 - 0.2e1 * t39 * t57;
t33 = -t60 * pkin(3) - pkin(4);
t31 = t37 * pkin(3) + qJ(5);
t13 = -pkin(4) * t45 - t24 * qJ(5) + t34;
t5 = t21 * pkin(4) - t22 * qJ(5) - t24 * qJD(5) + t35;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, t42; 0, 0, -t36, -t57, 0, 0, 0, 0, 0, -t40 * t36 - t54, t38 * t36 - t41 * t58, t46, -pkin(3) * t54 + t34 * t36 + t47, -t41 * t21 - t36 * t45, t46, t41 * t22 - t24 * t36, t13 * t36 - t41 * t5 + t47; 0, 0, 0, 0, 0.2e1 * t38 * t58, 0.2e1 * (-t38 ^ 2 + t40 ^ 2) * qJD(3), 0, 0, 0, t38 * t56, t40 * t56, t43, 0.2e1 * t34 * t35 + 0.2e1 * t52, 0.2e1 * t13 * t21 - 0.2e1 * t45 * t5, t43, -0.2e1 * t13 * t22 - 0.2e1 * t5 * t24, 0.2e1 * t13 * t5 + 0.2e1 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38 * t57 - t39 * t58, t64, 0, (t10 * t37 + t60 * t9) * pkin(3), t9, 0, t10, t19 * qJD(5) + t10 * t31 - t9 * t33; 0, 0, 0, 0, 0, 0, t58, -t59, 0, -pkin(6) * t58, pkin(6) * t59, (-t21 * t37 - t60 * t22) * pkin(3), (-t60 * t11 + t12 * t37) * pkin(3), -t11, qJD(5) * t45 - t31 * t21 + t33 * t22, t12, t16 * qJD(5) + t11 * t33 + t12 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t31 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t21, 0, -t22, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
