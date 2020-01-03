% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:33
% EndTime: 2019-12-31 17:03:33
% DurationCPUTime: 0.12s
% Computational Cost: add. (49->30), mult. (138->47), div. (0->0), fcn. (68->4), ass. (0->23)
t18 = 2 * qJD(3);
t13 = sin(qJ(4));
t15 = cos(qJ(4));
t22 = t15 * qJD(4);
t16 = cos(qJ(2));
t24 = pkin(1) * qJD(2);
t20 = t16 * t24;
t6 = qJD(3) + t20;
t14 = sin(qJ(2));
t8 = t14 * pkin(1) + qJ(3);
t26 = t6 * t13 + t8 * t22;
t21 = qJ(3) * qJD(4);
t25 = qJD(3) * t13 + t15 * t21;
t23 = t13 * qJD(4);
t10 = t14 * t24;
t19 = -t16 * pkin(1) - pkin(2);
t17 = -pkin(2) - pkin(6);
t12 = qJD(3) * t15;
t7 = -pkin(6) + t19;
t5 = -0.2e1 * t13 * t22;
t3 = t6 * t15;
t1 = 0.2e1 * (t13 ^ 2 - t15 ^ 2) * qJD(4);
t2 = [0, 0, 0, 0, -0.2e1 * t10, -0.2e1 * t20, 0.2e1 * t10, 0.2e1 * t6, 0.2e1 * t19 * t10 + 0.2e1 * t8 * t6, t5, t1, 0, 0, 0, 0.2e1 * t26, -0.2e1 * t8 * t23 + 0.2e1 * t3; 0, 0, 0, 0, -t10, -t20, t10, t18 + t20, -pkin(2) * t10 + t6 * qJ(3) + t8 * qJD(3), t5, t1, 0, 0, 0, t25 + t26, t12 + t3 + (-qJ(3) - t8) * t23; 0, 0, 0, 0, 0, 0, 0, t18, qJ(3) * t18, t5, t1, 0, 0, 0, 0.2e1 * t25, -0.2e1 * t13 * t21 + 0.2e1 * t12; 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t22, 0, t15 * t10 - t7 * t23, -t13 * t10 - t7 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t22, 0, -t17 * t23, -t17 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
