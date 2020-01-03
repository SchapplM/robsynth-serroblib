% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:50
% EndTime: 2019-12-31 18:14:50
% DurationCPUTime: 0.27s
% Computational Cost: add. (357->56), mult. (705->92), div. (0->0), fcn. (579->4), ass. (0->39)
t33 = sin(qJ(3));
t34 = cos(qJ(3));
t51 = sin(pkin(7));
t43 = qJD(3) * t51;
t52 = cos(pkin(7));
t44 = qJD(3) * t52;
t16 = t33 * t43 - t34 * t44;
t17 = -t33 * t44 - t34 * t43;
t46 = t52 * t34;
t19 = -t51 * t33 + t46;
t45 = t51 * t34;
t20 = -t52 * t33 - t45;
t63 = -0.2e1 * t16 * t20 - 0.2e1 * t19 * t17;
t62 = (-t51 * t16 + t52 * t17) * pkin(3);
t35 = -pkin(1) - pkin(6);
t53 = qJ(4) - t35;
t21 = t53 * t33;
t11 = -t51 * t21 + t53 * t46;
t12 = -t52 * t21 - t53 * t45;
t49 = t34 * qJD(3);
t15 = -t33 * qJD(4) - t53 * t49;
t50 = t33 * qJD(3);
t38 = -t34 * qJD(4) + t53 * t50;
t8 = t51 * t15 - t52 * t38;
t9 = t52 * t15 + t51 * t38;
t60 = t11 * t17 + t12 * t16 + t8 * t19 + t9 * t20;
t59 = 2 * qJD(2);
t58 = 2 * qJD(5);
t54 = t33 * pkin(3) + qJ(2);
t22 = pkin(3) * t49 + qJD(2);
t48 = qJ(2) * qJD(3);
t47 = t11 * t8 + t12 * t9;
t26 = t51 * pkin(3) + qJ(5);
t30 = -t52 * pkin(3) - pkin(4);
t39 = qJD(5) * t20 + t26 * t16 + t30 * t17;
t36 = 0.2e1 * t60;
t10 = -t20 * pkin(4) - t19 * qJ(5) + t54;
t5 = -t16 * pkin(4) - t17 * qJ(5) - t19 * qJD(5) + t22;
t1 = [0, 0, 0, 0, t59, qJ(2) * t59, -0.2e1 * t33 * t49, 0.2e1 * (t33 ^ 2 - t34 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t33 + 0.2e1 * t34 * t48, 0.2e1 * qJD(2) * t34 - 0.2e1 * t33 * t48, t36, 0.2e1 * t54 * t22 + 0.2e1 * t47, -0.2e1 * t10 * t16 - 0.2e1 * t5 * t20, t36, -0.2e1 * t10 * t17 - 0.2e1 * t5 * t19, 0.2e1 * t10 * t5 + 0.2e1 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t60, 0, t63, 0, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, 0, 0, 0, -t63; 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t49, 0, -t35 * t50, -t35 * t49, -t62, (t51 * t9 - t52 * t8) * pkin(3), -t8, t39, t9, t12 * qJD(5) + t9 * t26 + t8 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t49, 0, t62, t17, 0, -t16, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t26 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t16, 0, -t17, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
