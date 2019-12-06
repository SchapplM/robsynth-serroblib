% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:10
% EndTime: 2019-12-05 15:09:11
% DurationCPUTime: 0.15s
% Computational Cost: add. (88->34), mult. (271->63), div. (0->0), fcn. (229->6), ass. (0->29)
t16 = sin(pkin(8));
t18 = sin(qJ(3));
t29 = cos(pkin(8));
t31 = cos(qJ(3));
t8 = t18 * t16 - t31 * t29;
t32 = 2 * qJD(5);
t17 = sin(qJ(4));
t28 = qJD(4) * t17;
t19 = cos(qJ(4));
t13 = qJD(4) * t19;
t27 = qJD(4) * qJ(5);
t26 = -0.2e1 * pkin(3) * qJD(4);
t25 = pkin(6) * t28;
t24 = pkin(6) * t13;
t14 = t17 ^ 2;
t15 = t19 ^ 2;
t6 = t8 * qJD(3);
t23 = (-t14 - t15) * t6;
t21 = -t19 * pkin(4) - t17 * qJ(5);
t9 = t31 * t16 + t18 * t29;
t20 = t21 * qJD(4) + t19 * qJD(5);
t10 = -pkin(3) + t21;
t7 = t9 * qJD(3);
t5 = -pkin(4) * t28 + t17 * qJD(5) + t19 * t27;
t4 = -t7 * t19 + t8 * t28;
t3 = t8 * t13 + t7 * t17;
t2 = t9 * t13 - t17 * t6;
t1 = t19 * t6 + t9 * t28;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t9 * t23 + 0.2e1 * t8 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t7, t6, 0, 0, 0, 0, 0, t4, t3, t4, t23, -t3, pkin(6) * t23 + t7 * t10 - t8 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0.2e1 * t17 * t13, 0.2e1 * (-t14 + t15) * qJD(4), 0, 0, 0, t17 * t26, t19 * t26, 0.2e1 * t10 * t28 + 0.2e1 * t5 * t19, 0, -0.2e1 * t10 * t13 + 0.2e1 * t5 * t17, -0.2e1 * t10 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t2, 0, -t1, (pkin(4) * t6 - t9 * t27) * t17 + (-qJ(5) * t6 + (-qJD(4) * pkin(4) + qJD(5)) * t9) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t13, -t28, 0, t13, t5; 0, 0, 0, 0, 0, 0, 0, t13, -t28, 0, -t24, t25, -t24, t20, -t25, t20 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, qJ(5) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
