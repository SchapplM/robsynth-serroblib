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
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:33
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
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
% StartTime: 2021-01-15 15:32:32
% EndTime: 2021-01-15 15:32:35
% DurationCPUTime: 0.33s
% Computational Cost: add. (379->71), mult. (1049->136), div. (0->0), fcn. (913->6), ass. (0->51)
t41 = sin(qJ(2));
t42 = cos(qJ(3));
t43 = cos(qJ(2));
t62 = t43 * qJD(2);
t40 = sin(qJ(3));
t64 = t40 * qJD(3);
t69 = t41 * t64 - t42 * t62;
t68 = 2 * qJD(5);
t39 = sin(pkin(8));
t67 = t39 * t40;
t66 = qJ(4) + pkin(6);
t65 = cos(pkin(8));
t38 = t41 * qJD(2);
t63 = t42 * qJD(3);
t61 = -0.2e1 * pkin(2) * qJD(3);
t37 = pkin(3) * t64;
t59 = t43 * t64;
t58 = t40 * t62;
t36 = -t42 * pkin(3) - pkin(2);
t53 = qJD(3) * t66;
t21 = t42 * qJD(4) - t40 * t53;
t46 = -t40 * qJD(4) - t42 * t53;
t11 = t39 * t21 - t65 * t46;
t12 = t65 * t21 + t39 * t46;
t31 = t66 * t42;
t55 = t65 * t40;
t16 = t39 * t31 + t66 * t55;
t17 = t65 * t31 - t66 * t67;
t56 = t16 * t11 + t17 * t12;
t54 = t65 * t42;
t52 = t65 * t62;
t51 = qJD(3) * t54;
t23 = -t39 * t64 + t51;
t25 = t39 * t42 + t55;
t50 = -t43 * t23 + t25 * t38;
t22 = t25 * qJD(3);
t10 = t41 * t22 + t39 * t58 - t42 * t52;
t19 = t25 * t41;
t47 = t54 - t67;
t20 = t47 * t41;
t9 = t69 * t39 - t40 * t52 - t41 * t51;
t49 = -t10 * t17 + t19 * t11 + t20 * t12 - t9 * t16;
t48 = -t10 * t47 + t19 * t23 - t20 * t22 - t9 * t25;
t45 = 0.2e1 * t11 * t25 + 0.2e1 * t12 * t47 + 0.2e1 * t16 * t23 - 0.2e1 * t17 * t22;
t44 = -0.2e1 * t20 * t10 - 0.2e1 * t19 * t9 - 0.2e1 * t41 * t62;
t35 = -t65 * pkin(3) - pkin(4);
t33 = t39 * pkin(3) + qJ(5);
t14 = -pkin(4) * t47 - t25 * qJ(5) + t36;
t13 = -t43 * t22 - t38 * t47;
t5 = t22 * pkin(4) - t23 * qJ(5) - t25 * qJD(5) + t37;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, t44; 0, 0, -t38, -t62, 0, 0, 0, 0, 0, -t42 * t38 - t59, t40 * t38 - t43 * t63, t13, t50, t48, -pkin(3) * t59 + t36 * t38 + t49, t13, t48, -t50, t14 * t38 - t43 * t5 + t49; 0, 0, 0, 0, 0.2e1 * t40 * t63, 0.2e1 * (-t40 ^ 2 + t42 ^ 2) * qJD(3), 0, 0, 0, t40 * t61, t42 * t61, 0.2e1 * t36 * t22 - 0.2e1 * t37 * t47, 0.2e1 * t36 * t23 + 0.2e1 * t25 * t37, t45, 0.2e1 * t36 * t37 + 0.2e1 * t56, 0.2e1 * t14 * t22 - 0.2e1 * t47 * t5, t45, -0.2e1 * t14 * t23 - 0.2e1 * t5 * t25, 0.2e1 * t14 * t5 + 0.2e1 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41 * t63 - t58, t69, t9, t10, 0, (-t10 * t39 + t65 * t9) * pkin(3), t9, 0, -t10, t20 * qJD(5) - t10 * t33 - t9 * t35; 0, 0, 0, 0, 0, 0, t63, -t64, 0, -pkin(6) * t63, pkin(6) * t64, -t11, -t12, (-t22 * t39 - t65 * t23) * pkin(3), (-t65 * t11 + t12 * t39) * pkin(3), -t11, qJD(5) * t47 - t33 * t22 + t35 * t23, t12, t17 * qJD(5) + t11 * t35 + t12 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t33 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t23, 0, t37, t22, 0, -t23, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
