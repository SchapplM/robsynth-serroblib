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
% MMD_reg [((5+1)*5/2)x14]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 15:33:41
% EndTime: 2019-12-05 15:33:41
% DurationCPUTime: 0.17s
% Computational Cost: add. (122->41), mult. (307->84), div. (0->0), fcn. (255->6), ass. (0->29)
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t17 = sin(qJ(2));
t19 = cos(qJ(2));
t8 = t14 * t19 + t15 * t17;
t3 = t8 * qJD(2);
t7 = t14 * t17 - t15 * t19;
t28 = t7 * t3;
t10 = t14 * pkin(2) + pkin(6);
t27 = qJ(5) + t10;
t16 = sin(qJ(4));
t26 = t16 * qJD(4);
t18 = cos(qJ(4));
t25 = t18 * qJD(4);
t24 = 0.2e1 * t25;
t23 = pkin(4) * t26;
t11 = -t15 * pkin(2) - pkin(3);
t12 = t16 ^ 2;
t13 = t18 ^ 2;
t4 = t7 * qJD(2);
t22 = (t12 + t13) * t4;
t21 = qJD(4) * t27;
t20 = t16 * t4 - t8 * t25;
t9 = -t18 * pkin(4) + t11;
t6 = t27 * t18;
t5 = t27 * t16;
t2 = -t16 * qJD(5) - t18 * t21;
t1 = t18 * qJD(5) - t16 * t21;
t29 = [0, 0, 0, 0, -0.2e1 * t8 * t4 + 0.2e1 * t28, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t8 * t22 + 0.2e1 * t28; 0, 0, -t17 * qJD(2), -t19 * qJD(2), (-t14 * t4 - t15 * t3) * pkin(2), 0, 0, 0, 0, 0, -t3 * t18 + t7 * t26, t3 * t16 + t7 * t25, -t22, t3 * t9 + (-t4 * t6 + (qJD(4) * t5 + t1) * t8) * t18 + (-t2 * t8 - t4 * t5 + (pkin(4) * t7 - t6 * t8) * qJD(4)) * t16; 0, 0, 0, 0, 0, t16 * t24, 0.2e1 * (-t12 + t13) * qJD(4), 0, 0, 0, 0.2e1 * t11 * t26, t11 * t24, 0.2e1 * t1 * t18 - 0.2e1 * t2 * t16 + 0.2e1 * (-t16 * t6 + t18 * t5) * qJD(4), 0.2e1 * t6 * t1 - 0.2e1 * t5 * t2 + 0.2e1 * t9 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t16 + t2 * t18 + (t16 * t5 + t18 * t6) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t18 * t4 + t8 * t26, 0, t20 * pkin(4); 0, 0, 0, 0, 0, 0, 0, t25, -t26, 0, -t10 * t25, t10 * t26, -pkin(4) * t25, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t25, 0, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t29;
