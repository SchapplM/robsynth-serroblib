% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPR8
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
% MMD_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:23
% EndTime: 2019-12-31 17:08:24
% DurationCPUTime: 0.19s
% Computational Cost: add. (112->45), mult. (295->86), div. (0->0), fcn. (209->4), ass. (0->36)
t38 = 2 * qJD(3);
t37 = pkin(2) + pkin(3);
t36 = pkin(5) - pkin(6);
t24 = cos(qJ(2));
t35 = qJ(3) * t24;
t22 = sin(qJ(2));
t34 = t22 * qJ(3);
t33 = qJD(2) * t22;
t20 = qJD(2) * t24;
t21 = sin(qJ(4));
t32 = qJD(4) * t21;
t23 = cos(qJ(4));
t31 = qJD(4) * t23;
t30 = t22 * qJD(3);
t29 = -0.2e1 * pkin(1) * qJD(2);
t28 = pkin(5) * t33;
t27 = -t24 * pkin(2) - t34;
t10 = t22 * t21 + t24 * t23;
t26 = t27 * qJD(2) + t24 * qJD(3);
t19 = pkin(5) * t20;
t16 = t36 * t24;
t15 = t36 * t22;
t14 = -pkin(1) + t27;
t13 = -pkin(6) * t20 + t19;
t12 = t36 * t33;
t11 = -t24 * t21 + t22 * t23;
t9 = t37 * t24 + pkin(1) + t34;
t8 = -t30 + (pkin(2) * t22 - t35) * qJD(2);
t7 = t21 * qJD(3) + (qJ(3) * t23 - t21 * t37) * qJD(4);
t6 = -t23 * qJD(3) + (qJ(3) * t21 + t23 * t37) * qJD(4);
t5 = t30 + (-t37 * t22 + t35) * qJD(2);
t4 = (qJD(2) - qJD(4)) * t10;
t3 = t21 * t20 + t22 * t31 - t23 * t33 - t24 * t32;
t2 = -t21 * t12 - t23 * t13 + (t15 * t21 + t16 * t23) * qJD(4);
t1 = t23 * t12 - t21 * t13 + (-t15 * t23 + t16 * t21) * qJD(4);
t17 = [0, 0, 0, 0.2e1 * t22 * t20, 0.2e1 * (-t22 ^ 2 + t24 ^ 2) * qJD(2), 0, 0, 0, t22 * t29, t24 * t29, 0.2e1 * t14 * t33 - 0.2e1 * t8 * t24, 0, -0.2e1 * t14 * t20 - 0.2e1 * t8 * t22, 0.2e1 * t14 * t8, 0.2e1 * t11 * t4, -0.2e1 * t4 * t10 - 0.2e1 * t11 * t3, 0, 0, 0, 0.2e1 * t5 * t10 + 0.2e1 * t9 * t3, 0.2e1 * t5 * t11 + 0.2e1 * t9 * t4; 0, 0, 0, 0, 0, t20, -t33, 0, -t19, t28, -t19, t26, -t28, t26 * pkin(5), 0, 0, -t4, t3, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, qJ(3) * t38, 0, 0, 0, 0, 0, 0.2e1 * t7, -0.2e1 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, t19, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
