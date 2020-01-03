% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:37
% EndTime: 2019-12-31 17:00:37
% DurationCPUTime: 0.14s
% Computational Cost: add. (78->37), mult. (206->75), div. (0->0), fcn. (109->2), ass. (0->26)
t15 = cos(qJ(2));
t11 = t15 * qJD(3);
t14 = sin(qJ(2));
t24 = t14 * qJ(3);
t18 = -t15 * pkin(2) - t24;
t28 = t18 * qJD(2) + t11;
t13 = -pkin(2) - qJ(4);
t27 = t13 * t15 - t24;
t26 = pkin(3) + pkin(5);
t23 = qJD(2) * t14;
t12 = qJD(2) * t15;
t22 = qJ(3) * qJD(3);
t21 = -0.2e1 * pkin(1) * qJD(2);
t20 = pkin(5) * t23;
t19 = pkin(2) * t23 - t14 * qJD(3);
t16 = 0.2e1 * qJD(3);
t9 = pkin(5) * t12;
t8 = t26 * t15;
t7 = t26 * t14;
t6 = -pkin(1) + t18;
t5 = pkin(3) * t12 + t9;
t4 = t26 * t23;
t3 = -pkin(1) + t27;
t2 = -qJ(3) * t12 + t19;
t1 = -t15 * qJD(4) + (-qJ(3) * t15 + qJ(4) * t14) * qJD(2) + t19;
t10 = [0, 0, 0, 0.2e1 * t14 * t12, 0.2e1 * (-t14 ^ 2 + t15 ^ 2) * qJD(2), 0, 0, 0, t14 * t21, t15 * t21, 0, 0.2e1 * t2 * t15 - 0.2e1 * t6 * t23, -0.2e1 * t6 * t12 - 0.2e1 * t2 * t14, 0.2e1 * t6 * t2, 0.2e1 * t5 * t14 - 0.2e1 * t4 * t15 + 0.2e1 * (-t14 * t8 + t15 * t7) * qJD(2), -0.2e1 * t1 * t14 - 0.2e1 * t3 * t12, -0.2e1 * t1 * t15 + 0.2e1 * t3 * t23, 0.2e1 * t3 * t1 - 0.2e1 * t8 * t4 + 0.2e1 * t7 * t5; 0, 0, 0, 0, 0, t12, -t23, 0, -t9, t20, t28, t9, -t20, t28 * pkin(5), t27 * qJD(2) - qJD(4) * t14 + t11, -t4, -t5, -t4 * qJ(3) + t8 * qJD(3) - t7 * qJD(4) + t5 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0.2e1 * t22, 0, t16, 0.2e1 * qJD(4), -0.2e1 * t13 * qJD(4) + 0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, t9, t12, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
