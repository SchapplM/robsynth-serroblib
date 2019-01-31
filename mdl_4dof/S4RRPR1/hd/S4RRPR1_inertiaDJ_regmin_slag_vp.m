% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x10]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-01-31 13:16
% Revision: 9ef80adae39e3cd5824e7abdb6e4e1e7895c437e (2019-01-31)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-31 13:16:45
% EndTime: 2019-01-31 13:16:46
% DurationCPUTime: 0.11s
% Computational Cost: add. (86->26), mult. (254->50), div. (0->0), fcn. (174->6), ass. (0->29)
t15 = sin(pkin(7));
t31 = pkin(2) * t15;
t16 = cos(pkin(7));
t13 = t16 * pkin(2) + pkin(3);
t20 = cos(qJ(2));
t14 = t20 * pkin(1) + pkin(2);
t18 = sin(qJ(2));
t29 = t15 * t18;
t21 = -pkin(1) * t29 + t16 * t14;
t7 = pkin(3) + t21;
t30 = -t13 - t7;
t28 = t16 * t18;
t27 = pkin(1) * qJD(2);
t17 = sin(qJ(4));
t26 = qJD(4) * t17;
t25 = t18 * t27;
t24 = t20 * t27;
t10 = pkin(1) * t28 + t15 * t14;
t8 = (-t15 * t20 - t28) * t27;
t23 = t10 * t26 - t17 * t8;
t19 = cos(qJ(4));
t9 = (t16 * t20 - t29) * t27;
t22 = -t17 * t9 + t19 * t8;
t12 = t26 * t31;
t6 = (-t13 * t17 - t19 * t31) * qJD(4);
t5 = -qJD(4) * t19 * t13 + t12;
t2 = (-t10 * t19 - t17 * t7) * qJD(4) + t22;
t1 = (-qJD(4) * t7 - t9) * t19 + t23;
t3 = [0, 0, 0, 0, -0.2e1 * t25, -0.2e1 * t24, 0.2e1 * t10 * t9 + 0.2e1 * t21 * t8, 0, 0.2e1 * t2, 0.2e1 * t1; 0, 0, 0, 0, -t25, -t24 (t15 * t9 + t16 * t8) * pkin(2), 0 ((-t10 - t31) * t19 + t30 * t17) * qJD(4) + t22, t12 + (t30 * qJD(4) - t9) * t19 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t6, 0.2e1 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
