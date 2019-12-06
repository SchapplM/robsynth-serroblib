% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:15
% EndTime: 2019-12-05 17:38:16
% DurationCPUTime: 0.19s
% Computational Cost: add. (130->44), mult. (280->70), div. (0->0), fcn. (222->4), ass. (0->31)
t32 = qJD(4) + qJD(5);
t16 = -pkin(6) + qJ(2);
t31 = pkin(7) - t16;
t20 = cos(qJ(5));
t21 = cos(qJ(4));
t30 = t20 * t21;
t17 = pkin(1) + qJ(3);
t18 = sin(qJ(5));
t29 = qJD(5) * t18;
t19 = sin(qJ(4));
t28 = t19 * qJD(2);
t27 = t19 * qJD(4);
t26 = t21 * qJD(4);
t25 = qJ(2) * qJD(2);
t24 = pkin(4) * t29;
t23 = qJD(5) * t20 * pkin(4);
t10 = t31 * t21;
t7 = t18 * t21 + t20 * t19;
t22 = 0.2e1 * qJD(2);
t15 = t21 * qJD(2);
t14 = t19 * pkin(4) + t17;
t11 = pkin(4) * t26 + qJD(3);
t9 = t31 * t19;
t8 = -t18 * t19 + t30;
t6 = -qJD(4) * t10 + t28;
t5 = t31 * t27 + t15;
t4 = -t18 * t27 - t19 * t29 + t32 * t30;
t3 = t32 * t7;
t2 = -t18 * t6 + t20 * t5 + (t10 * t18 + t20 * t9) * qJD(5);
t1 = -t18 * t5 - t20 * t6 + (t10 * t20 - t18 * t9) * qJD(5);
t12 = [0, 0, 0, 0, t22, 0.2e1 * t25, t22, 0.2e1 * qJD(3), 0.2e1 * t17 * qJD(3) + 0.2e1 * t25, -0.2e1 * t19 * t26, 0.2e1 * (t19 ^ 2 - t21 ^ 2) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t19 + 0.2e1 * t17 * t26, 0.2e1 * qJD(3) * t21 - 0.2e1 * t17 * t27, -0.2e1 * t8 * t3, 0.2e1 * t3 * t7 - 0.2e1 * t8 * t4, 0, 0, 0, 0.2e1 * t11 * t7 + 0.2e1 * t14 * t4, 0.2e1 * t11 * t8 - 0.2e1 * t14 * t3; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, 0, 0, 0, 0, -t26, t27, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t26, 0, -t16 * t27 + t15, -t16 * t26 - t28, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t26, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t24, -0.2e1 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
