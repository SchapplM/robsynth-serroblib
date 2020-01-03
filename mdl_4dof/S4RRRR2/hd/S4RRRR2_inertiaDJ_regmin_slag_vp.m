% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:18
% EndTime: 2019-12-31 17:23:19
% DurationCPUTime: 0.26s
% Computational Cost: add. (225->55), mult. (601->91), div. (0->0), fcn. (474->6), ass. (0->58)
t68 = qJD(3) + qJD(4);
t67 = -pkin(6) - pkin(7);
t44 = cos(qJ(2));
t66 = t44 * pkin(1);
t41 = sin(qJ(2));
t34 = t41 * pkin(1) + pkin(6);
t65 = -pkin(7) - t34;
t39 = sin(qJ(4));
t40 = sin(qJ(3));
t42 = cos(qJ(4));
t43 = cos(qJ(3));
t23 = t39 * t43 + t42 * t40;
t12 = t68 * t23;
t22 = t39 * t40 - t42 * t43;
t57 = t40 * qJD(3);
t52 = pkin(3) * t57;
t59 = pkin(1) * qJD(2);
t54 = t41 * t59;
t26 = t52 + t54;
t36 = -t43 * pkin(3) - pkin(2);
t29 = t36 - t66;
t64 = t29 * t12 + t26 * t22;
t11 = t68 * t22;
t63 = -t29 * t11 + t26 * t23;
t62 = -t36 * t11 + t23 * t52;
t61 = t36 * t12 + t22 * t52;
t35 = -pkin(2) - t66;
t37 = t43 * qJD(3);
t60 = t35 * t37 + t40 * t54;
t58 = pkin(3) * qJD(4);
t56 = pkin(2) * t57;
t55 = pkin(2) * t37;
t53 = t44 * t59;
t51 = t39 * t58;
t50 = t42 * t58;
t49 = qJD(3) * t67;
t48 = qJD(3) * t65;
t47 = t40 * t53;
t46 = t43 * t53;
t45 = t35 * t57 - t43 * t54;
t38 = t43 * pkin(7);
t33 = 0.2e1 * t40 * t37;
t31 = t43 * pkin(6) + t38;
t30 = t67 * t40;
t25 = t43 * t49;
t24 = t40 * t49;
t21 = 0.2e1 * (-t40 ^ 2 + t43 ^ 2) * qJD(3);
t20 = t43 * t34 + t38;
t19 = t65 * t40;
t16 = t43 * t48 - t47;
t15 = t40 * t48 + t46;
t6 = -0.2e1 * t23 * t11;
t5 = -t39 * t24 + t42 * t25 + (-t30 * t39 - t31 * t42) * qJD(4);
t4 = -t42 * t24 - t39 * t25 + (-t30 * t42 + t31 * t39) * qJD(4);
t3 = 0.2e1 * t11 * t22 - 0.2e1 * t23 * t12;
t2 = -t39 * t15 + t42 * t16 + (-t19 * t39 - t20 * t42) * qJD(4);
t1 = -t42 * t15 - t39 * t16 + (-t19 * t42 + t20 * t39) * qJD(4);
t7 = [0, 0, 0, 0, -0.2e1 * t54, -0.2e1 * t53, t33, t21, 0, 0, 0, 0.2e1 * t45, 0.2e1 * t60, t6, t3, 0, 0, 0, 0.2e1 * t64, 0.2e1 * t63; 0, 0, 0, 0, -t54, -t53, t33, t21, 0, 0, 0, t45 - t56, -t55 + t60, t6, t3, 0, 0, 0, t61 + t64, t62 + t63; 0, 0, 0, 0, 0, 0, t33, t21, 0, 0, 0, -0.2e1 * t56, -0.2e1 * t55, t6, t3, 0, 0, 0, 0.2e1 * t61, 0.2e1 * t62; 0, 0, 0, 0, 0, 0, 0, 0, t37, -t57, 0, -t34 * t37 - t47, t34 * t57 - t46, 0, 0, -t11, -t12, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, t37, -t57, 0, -pkin(6) * t37, pkin(6) * t57, 0, 0, -t11, -t12, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t51, -0.2e1 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
