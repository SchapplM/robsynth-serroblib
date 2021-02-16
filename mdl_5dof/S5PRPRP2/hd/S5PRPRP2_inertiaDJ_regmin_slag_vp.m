% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRP2
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
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:52
% EndTime: 2021-01-15 15:04:53
% DurationCPUTime: 0.25s
% Computational Cost: add. (172->53), mult. (458->99), div. (0->0), fcn. (346->4), ass. (0->38)
t16 = sin(pkin(8));
t43 = (-qJ(5) - pkin(6)) * t16;
t17 = cos(pkin(8));
t18 = sin(qJ(4));
t19 = cos(qJ(4));
t33 = t16 * qJD(5);
t37 = t18 * qJ(3);
t40 = t17 * pkin(3);
t27 = pkin(2) + t40;
t20 = -t16 * pkin(6) - t27;
t34 = qJD(4) * t19;
t36 = qJD(3) * t19;
t39 = -t17 * t36 - t20 * t34;
t1 = t18 * t33 + (t19 * t16 * qJ(5) + t17 * t37) * qJD(4) + t39;
t32 = t18 * qJD(3);
t28 = t17 * t32;
t30 = t19 * t17 * qJ(3);
t2 = -t28 - t19 * t33 + (-t30 + (t27 - t43) * t18) * qJD(4);
t21 = -pkin(2) + t43;
t3 = t21 * t19 + (-t19 * pkin(3) - pkin(4) - t37) * t17;
t6 = t30 + (t21 - t40) * t18;
t42 = t1 * t18 - t2 * t19 + (t18 * t3 - t19 * t6) * qJD(4);
t41 = 0.2e1 * qJD(4);
t35 = qJD(4) * t18;
t31 = qJ(3) * qJD(4);
t11 = t16 * t35;
t29 = t16 * t34;
t26 = t18 * t31;
t25 = t16 * t17 * t41;
t14 = t16 ^ 2;
t24 = 0.2e1 * (t17 ^ 2 + t14) * qJD(3);
t13 = t17 * t34;
t12 = t17 * t35;
t9 = (pkin(4) * t18 + qJ(3)) * t16;
t8 = (pkin(4) * t34 + qJD(3)) * t16;
t5 = -t28 + (-t18 * t20 - t30) * qJD(4);
t4 = t17 * t26 + t39;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * t8 + (-t1 * t19 - t18 * t2 + (-t18 * t6 - t19 * t3) * qJD(4)) * t16; 0, 0, 0, 0, 0, t24, qJ(3) * t24, -0.2e1 * t14 * t18 * t34, (t18 ^ 2 - t19 ^ 2) * t14 * t41, t18 * t25, t19 * t25, 0, -0.2e1 * t5 * t17 + 0.2e1 * (t19 * t31 + t32) * t14, -0.2e1 * t4 * t17 + 0.2e1 * (-t26 + t36) * t14, -0.2e1 * t2 * t17 + 0.2e1 * (t18 * t8 + t9 * t34) * t16, -0.2e1 * t1 * t17 + 0.2e1 * (t19 * t8 - t9 * t35) * t16, 0.2e1 * t42 * t16, -0.2e1 * t6 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t9 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t13, t12, t13, 0, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t11, -t29, t11, 0, -pkin(4) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t29, 0, t5, t4, t2, t1, pkin(4) * t11, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t34, -t35, -t34, 0, -pkin(4) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t11, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
