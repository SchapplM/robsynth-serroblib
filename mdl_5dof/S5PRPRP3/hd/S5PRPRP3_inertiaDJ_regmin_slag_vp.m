% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:14:01
% EndTime: 2021-01-15 15:14:02
% DurationCPUTime: 0.17s
% Computational Cost: add. (148->42), mult. (377->81), div. (0->0), fcn. (316->6), ass. (0->36)
t36 = 2 * qJD(4);
t20 = sin(pkin(8));
t21 = cos(pkin(8));
t23 = sin(qJ(2));
t25 = cos(qJ(2));
t11 = t20 * t23 - t21 * t25;
t12 = t20 * t25 + t21 * t23;
t7 = t12 * qJD(2);
t35 = t11 * t7;
t24 = cos(qJ(4));
t34 = t24 * pkin(4);
t14 = t20 * pkin(2) + pkin(6);
t32 = qJ(5) + t14;
t10 = t32 * t24;
t33 = t10 * t24;
t22 = sin(qJ(4));
t31 = t22 * qJD(4);
t17 = t24 * qJD(4);
t30 = 0.2e1 * t31;
t29 = 0.2e1 * t17;
t28 = pkin(4) * t31;
t15 = -t21 * pkin(2) - pkin(3);
t18 = t22 ^ 2;
t19 = t24 ^ 2;
t8 = t11 * qJD(2);
t27 = (t18 + t19) * t8;
t2 = -t12 * t17 + t22 * t8;
t5 = -t24 * qJD(5) + t32 * t31;
t6 = -t22 * qJD(5) - t32 * t17;
t9 = t32 * t22;
t26 = -t6 * t22 - t5 * t24 + (-t10 * t22 + t24 * t9) * qJD(4);
t13 = t15 - t34;
t4 = t11 * t31 - t7 * t24;
t3 = t11 * t17 + t7 * t22;
t1 = t12 * t31 + t24 * t8;
t16 = [0, 0, 0, 0, -0.2e1 * t12 * t8 + 0.2e1 * t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t12 * t27 + 0.2e1 * t35; 0, 0, -t23 * qJD(2), -t25 * qJD(2), (-t20 * t8 - t21 * t7) * pkin(2), 0, 0, 0, 0, 0, t4, t3, t4, t3, -t27, -t8 * t33 + t7 * t13 + (pkin(4) * qJD(4) * t11 - t8 * t9) * t22 + t26 * t12; 0, 0, 0, 0, 0, t22 * t29, (-t18 + t19) * t36, 0, 0, 0, t15 * t30, t15 * t29, (t13 - t34) * t30, (pkin(4) * t18 + t13 * t24) * t36, 0.2e1 * t26, -0.2e1 * t10 * t5 + 0.2e1 * t13 * t28 - 0.2e1 * t9 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t22 + t6 * t24 + (t22 * t9 + t33) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, t17, -t31, 0, -t14 * t17, t14 * t31, t6, t5, -pkin(4) * t17, t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t17, -t31, -t17, 0, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t17, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t16;
