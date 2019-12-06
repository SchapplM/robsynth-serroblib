% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:45
% EndTime: 2019-12-05 16:00:46
% DurationCPUTime: 0.22s
% Computational Cost: add. (146->49), mult. (370->87), div. (0->0), fcn. (318->6), ass. (0->35)
t36 = qJD(4) + qJD(5);
t23 = cos(qJ(2));
t35 = t36 * t23;
t34 = 2 * qJD(3);
t24 = -pkin(2) - pkin(6);
t33 = pkin(7) - t24;
t21 = cos(qJ(5));
t22 = cos(qJ(4));
t32 = t21 * t22;
t18 = sin(qJ(5));
t31 = qJD(5) * t18;
t19 = sin(qJ(4));
t30 = t19 * qJD(4);
t20 = sin(qJ(2));
t17 = t20 * qJD(2);
t29 = t22 * qJD(4);
t28 = t23 * qJD(2);
t27 = qJ(3) * qJD(4);
t26 = pkin(4) * t31;
t25 = qJD(5) * t21 * pkin(4);
t12 = t33 * t22;
t9 = t18 * t22 + t21 * t19;
t10 = -t18 * t19 + t32;
t16 = t19 * pkin(4) + qJ(3);
t13 = pkin(4) * t29 + qJD(3);
t11 = t33 * t19;
t8 = qJD(4) * t12;
t7 = t33 * t30;
t6 = -t18 * t30 - t19 * t31 + t36 * t32;
t5 = t36 * t9;
t4 = t10 * t17 + t9 * t35;
t3 = t10 * t35 - t9 * t17;
t2 = t18 * t8 + t21 * t7 + (t11 * t21 + t12 * t18) * qJD(5);
t1 = -t18 * t7 + t21 * t8 + (-t11 * t18 + t12 * t21) * qJD(5);
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t17, -t28, t17, t28, t20 * qJD(3) + (-pkin(2) * t20 + qJ(3) * t23) * qJD(2), 0, 0, 0, 0, 0, t19 * t28 + t20 * t29, -t20 * t30 + t22 * t28, 0, 0, 0, 0, 0, t20 * t6 + t9 * t28, t10 * t28 - t20 * t5; 0, 0, 0, 0, 0, t34, qJ(3) * t34, -0.2e1 * t19 * t29, 0.2e1 * (t19 ^ 2 - t22 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t19 + 0.2e1 * t22 * t27, 0.2e1 * qJD(3) * t22 - 0.2e1 * t19 * t27, -0.2e1 * t10 * t5, -0.2e1 * t10 * t6 + 0.2e1 * t5 * t9, 0, 0, 0, 0.2e1 * t13 * t9 + 0.2e1 * t16 * t6, 0.2e1 * t13 * t10 - 0.2e1 * t16 * t5; 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t17 + t23 * t30, -t19 * t17 + t23 * t29, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t29, 0, -t24 * t30, -t24 * t29, 0, 0, -t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t29, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t26, -0.2e1 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
