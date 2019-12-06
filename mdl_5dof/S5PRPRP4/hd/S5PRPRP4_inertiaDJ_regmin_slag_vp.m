% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRP4
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:05
% EndTime: 2019-12-05 15:36:06
% DurationCPUTime: 0.16s
% Computational Cost: add. (112->38), mult. (309->69), div. (0->0), fcn. (265->6), ass. (0->31)
t31 = 2 * qJD(5);
t17 = sin(pkin(8));
t18 = cos(pkin(8));
t20 = sin(qJ(2));
t22 = cos(qJ(2));
t10 = t17 * t22 + t18 * t20;
t6 = t10 * qJD(2);
t9 = t17 * t20 - t18 * t22;
t30 = t9 * t6;
t19 = sin(qJ(4));
t29 = qJD(4) * t19;
t21 = cos(qJ(4));
t14 = qJD(4) * t21;
t28 = 0.2e1 * t14;
t11 = t17 * pkin(2) + pkin(6);
t27 = t11 * t29;
t26 = t11 * t14;
t12 = -t18 * pkin(2) - pkin(3);
t15 = t19 ^ 2;
t16 = t21 ^ 2;
t7 = t9 * qJD(2);
t25 = (t15 + t16) * t7;
t24 = -t21 * pkin(4) - t19 * qJ(5);
t23 = t24 * qJD(4) + t21 * qJD(5);
t8 = t12 + t24;
t5 = -pkin(4) * t29 + qJ(5) * t14 + t19 * qJD(5);
t4 = -t6 * t21 + t9 * t29;
t3 = t9 * t14 + t6 * t19;
t2 = t10 * t14 - t19 * t7;
t1 = t10 * t29 + t21 * t7;
t13 = [0, 0, 0, 0, -0.2e1 * t10 * t7 + 0.2e1 * t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t10 * t25 + 0.2e1 * t30; 0, 0, -qJD(2) * t20, -qJD(2) * t22, (-t17 * t7 - t18 * t6) * pkin(2), 0, 0, 0, 0, 0, t4, t3, t4, -t25, -t3, -t11 * t25 - t9 * t5 + t6 * t8; 0, 0, 0, 0, 0, t19 * t28, 0.2e1 * (-t15 + t16) * qJD(4), 0, 0, 0, 0.2e1 * t12 * t29, t12 * t28, 0.2e1 * t5 * t21 + 0.2e1 * t8 * t29, 0, -0.2e1 * t8 * t14 + 0.2e1 * t5 * t19, -0.2e1 * t8 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t2, 0, -t1, -(-pkin(4) * t19 + qJ(5) * t21) * t7 + t23 * t10; 0, 0, 0, 0, 0, 0, 0, t14, -t29, 0, -t26, t27, -t26, t23, -t27, t23 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t14, -t29, 0, t14, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, qJ(5) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
