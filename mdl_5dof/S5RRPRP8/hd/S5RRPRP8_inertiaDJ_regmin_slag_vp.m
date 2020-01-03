% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:26
% EndTime: 2019-12-31 20:04:27
% DurationCPUTime: 0.35s
% Computational Cost: add. (384->79), mult. (853->141), div. (0->0), fcn. (672->4), ass. (0->44)
t43 = sin(qJ(2));
t45 = cos(qJ(2));
t60 = -t45 * pkin(2) - t43 * qJ(3);
t52 = t43 * qJD(2);
t57 = pkin(6) - pkin(7);
t22 = t57 * t52;
t38 = t45 * qJD(2);
t35 = pkin(6) * t38;
t23 = -pkin(7) * t38 + t35;
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t26 = t57 * t43;
t27 = t57 * t45;
t47 = t42 * t26 + t44 * t27;
t4 = t47 * qJD(4) - t42 * t22 - t44 * t23;
t59 = 2 * qJD(3);
t58 = -pkin(2) - pkin(3);
t56 = qJ(3) * t38 + t43 * qJD(3);
t55 = qJD(4) * t42;
t54 = qJD(4) * t44;
t53 = qJD(4) * t45;
t51 = -0.2e1 * pkin(1) * qJD(2);
t25 = -pkin(1) + t60;
t50 = pkin(6) * t52;
t49 = t44 * t58;
t17 = t45 * pkin(3) - t25;
t19 = t43 * t42 + t45 * t44;
t24 = t44 * qJ(3) + t42 * t58;
t3 = t44 * t22 - t42 * t23 - t26 * t54 + t27 * t55;
t10 = t58 * t52 + t56;
t46 = t60 * qJD(2) + t45 * qJD(3);
t21 = -t42 * qJ(3) - pkin(4) + t49;
t20 = -t45 * t42 + t43 * t44;
t13 = pkin(2) * t52 - t56;
t12 = t42 * qJD(3) + qJD(4) * t24;
t11 = qJ(3) * t55 - t44 * qJD(3) - qJD(4) * t49;
t9 = t19 * qJD(2) - t43 * t55 - t44 * t53;
t8 = t43 * t54 - t44 * t52 + (t38 - t53) * t42;
t7 = -t19 * qJ(5) + t47;
t6 = -t20 * qJ(5) + t44 * t26 - t42 * t27;
t5 = t8 * pkin(4) + t10;
t2 = -t9 * qJ(5) - t20 * qJD(5) - t4;
t1 = -t8 * qJ(5) - t19 * qJD(5) - t3;
t14 = [0, 0, 0, 0.2e1 * t43 * t38, 0.2e1 * (-t43 ^ 2 + t45 ^ 2) * qJD(2), 0, 0, 0, t43 * t51, t45 * t51, -0.2e1 * t13 * t45 + 0.2e1 * t25 * t52, 0, -0.2e1 * t13 * t43 - 0.2e1 * t25 * t38, 0.2e1 * t25 * t13, 0.2e1 * t20 * t9, -0.2e1 * t9 * t19 - 0.2e1 * t20 * t8, 0, 0, 0, 0.2e1 * t10 * t19 + 0.2e1 * t17 * t8, 0.2e1 * t10 * t20 + 0.2e1 * t17 * t9, -0.2e1 * t1 * t19 - 0.2e1 * t2 * t20 - 0.2e1 * t6 * t9 - 0.2e1 * t7 * t8, 0.2e1 * t7 * t1 + 0.2e1 * t6 * t2 + 0.2e1 * (t19 * pkin(4) + t17) * t5; 0, 0, 0, 0, 0, t38, -t52, 0, -t35, t50, -t35, t46, -t50, t46 * pkin(6), 0, 0, -t9, t8, 0, t4, -t3, t11 * t19 + t12 * t20 - t21 * t9 - t24 * t8, t1 * t24 - t7 * t11 - t6 * t12 + t2 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, qJ(3) * t59, 0, 0, 0, 0, 0, 0.2e1 * t12, -0.2e1 * t11, 0, -0.2e1 * t24 * t11 - 0.2e1 * t21 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t35, 0, 0, 0, 0, 0, 0, 0, -t42 * t8 - t44 * t9 + (-t19 * t44 + t20 * t42) * qJD(4), t1 * t42 + t2 * t44 + (-t42 * t6 + t44 * t7) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t54, 0, -t11 * t42 - t12 * t44 + (-t21 * t42 + t24 * t44) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, -t4, t3, -pkin(4) * t9, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, -t12 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t54, 0, -pkin(4) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
