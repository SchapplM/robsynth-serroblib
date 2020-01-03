% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:37
% EndTime: 2019-12-31 17:38:38
% DurationCPUTime: 0.16s
% Computational Cost: add. (64->33), mult. (168->69), div. (0->0), fcn. (148->6), ass. (0->23)
t28 = 2 * qJD(3);
t14 = sin(pkin(8));
t15 = cos(pkin(8));
t20 = -pkin(2) - pkin(3);
t8 = t15 * qJ(3) + t14 * t20;
t17 = sin(qJ(2));
t27 = qJD(2) * t17;
t19 = cos(qJ(2));
t26 = qJD(2) * t19;
t16 = sin(qJ(5));
t25 = qJD(5) * t16;
t18 = cos(qJ(5));
t24 = qJD(5) * t18;
t23 = t14 * qJD(3);
t22 = t15 * qJD(3);
t21 = t17 * t14 + t19 * t15;
t7 = -t14 * qJ(3) + t15 * t20;
t6 = -pkin(6) + t8;
t5 = pkin(4) - t7;
t4 = -t19 * t14 + t17 * t15;
t2 = t21 * qJD(2);
t1 = -t14 * t26 + t15 * t27;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t1 * t21 + 0.2e1 * t4 * t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t27, -t26, -t27, t26, t17 * qJD(3) + (-pkin(2) * t17 + qJ(3) * t19) * qJD(2), -t1, t2, t1 * t7 + t2 * t8 + (t14 * t21 + t15 * t4) * qJD(3), 0, 0, 0, 0, 0, -t1 * t18 - t21 * t25, t1 * t16 - t21 * t24; 0, 0, 0, 0, 0, t28, qJ(3) * t28, 0.2e1 * t23, 0.2e1 * t22, (-t14 * t7 + t15 * t8) * t28, 0.2e1 * t16 * t24, 0.2e1 * (-t16 ^ 2 + t18 ^ 2) * qJD(5), 0, 0, 0, 0.2e1 * t18 * t23 - 0.2e1 * t5 * t25, -0.2e1 * t16 * t23 - 0.2e1 * t5 * t24; 0, 0, 0, 0, 0, 0, t27, 0, 0, t1 * t15 + t2 * t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t25, t15 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t2 - t4 * t24, -t18 * t2 + t4 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t25, 0, -t16 * t22 - t6 * t24, -t18 * t22 + t6 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t24, t14 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
