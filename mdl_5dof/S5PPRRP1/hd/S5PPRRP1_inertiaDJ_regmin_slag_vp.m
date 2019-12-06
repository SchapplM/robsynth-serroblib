% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRRP1
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
% MMD_reg [((5+1)*5/2)x14]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:14
% EndTime: 2019-12-05 15:07:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (89->35), mult. (260->74), div. (0->0), fcn. (210->6), ass. (0->26)
t13 = sin(pkin(8));
t15 = sin(qJ(3));
t25 = cos(pkin(8));
t28 = cos(qJ(3));
t5 = t15 * t13 - t28 * t25;
t26 = -qJ(5) - pkin(6);
t14 = sin(qJ(4));
t24 = t14 * qJD(4);
t16 = cos(qJ(4));
t23 = t16 * qJD(4);
t22 = -0.2e1 * pkin(3) * qJD(4);
t21 = pkin(4) * t24;
t11 = t14 ^ 2;
t12 = t16 ^ 2;
t3 = t5 * qJD(3);
t20 = (-t11 - t12) * t3;
t19 = qJD(4) * t26;
t6 = t28 * t13 + t15 * t25;
t17 = t14 * t3 - t6 * t23;
t10 = -t16 * pkin(4) - pkin(3);
t8 = t26 * t16;
t7 = t26 * t14;
t4 = t6 * qJD(3);
t2 = -t14 * qJD(5) + t16 * t19;
t1 = t16 * qJD(5) + t14 * t19;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t6 * t20 + 0.2e1 * t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t4, t3, 0, 0, 0, 0, 0, -t4 * t16 + t5 * t24, t4 * t14 + t5 * t23, t20, t4 * t10 + (t3 * t8 + (-qJD(4) * t7 + t1) * t6) * t16 + (-t2 * t6 + t3 * t7 + (pkin(4) * t5 + t6 * t8) * qJD(4)) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t1 + t16 * t2 + (-t14 * t7 - t16 * t8) * qJD(4); 0, 0, 0, 0, 0, 0.2e1 * t14 * t23, 0.2e1 * (-t11 + t12) * qJD(4), 0, 0, 0, t14 * t22, t16 * t22, 0.2e1 * t1 * t16 - 0.2e1 * t2 * t14 + 0.2e1 * (t14 * t8 - t16 * t7) * qJD(4), -0.2e1 * t8 * t1 + 0.2e1 * t10 * t21 + 0.2e1 * t7 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t16 * t3 + t6 * t24, 0, t17 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, 0, -t21; 0, 0, 0, 0, 0, 0, 0, t23, -t24, 0, -pkin(6) * t23, pkin(6) * t24, -pkin(4) * t23, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
