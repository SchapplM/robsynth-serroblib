% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:16
% DurationCPUTime: 0.14s
% Computational Cost: add. (108->37), mult. (246->65), div. (0->0), fcn. (194->4), ass. (0->29)
t30 = qJD(3) + qJD(4);
t29 = 2 * qJD(2);
t19 = -pkin(1) - pkin(5);
t28 = pkin(6) - t19;
t17 = cos(qJ(4));
t18 = cos(qJ(3));
t27 = t17 * t18;
t16 = sin(qJ(3));
t26 = qJD(3) * t16;
t25 = qJD(3) * t18;
t24 = qJD(3) * t19;
t15 = sin(qJ(4));
t23 = qJD(4) * t15;
t22 = qJ(2) * qJD(3);
t21 = pkin(3) * t23;
t20 = qJD(4) * t17 * pkin(3);
t10 = t28 * t18;
t7 = t15 * t18 + t17 * t16;
t14 = t16 * pkin(3) + qJ(2);
t11 = pkin(3) * t25 + qJD(2);
t9 = t28 * t16;
t8 = -t15 * t16 + t27;
t6 = qJD(3) * t10;
t5 = t28 * t26;
t4 = -t15 * t26 - t16 * t23 + t30 * t27;
t3 = t30 * t7;
t2 = t15 * t6 + t17 * t5 + (t10 * t15 + t17 * t9) * qJD(4);
t1 = -t15 * t5 + t17 * t6 + (t10 * t17 - t15 * t9) * qJD(4);
t12 = [0, 0, 0, 0, t29, qJ(2) * t29, -0.2e1 * t16 * t25, 0.2e1 * (t16 ^ 2 - t18 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t16 + 0.2e1 * t18 * t22, 0.2e1 * qJD(2) * t18 - 0.2e1 * t16 * t22, -0.2e1 * t8 * t3, 0.2e1 * t3 * t7 - 0.2e1 * t8 * t4, 0, 0, 0, 0.2e1 * t11 * t7 + 0.2e1 * t14 * t4, 0.2e1 * t11 * t8 - 0.2e1 * t14 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t25, 0, -t16 * t24, -t18 * t24, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t25, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t21, -0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
