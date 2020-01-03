% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:34
% DurationCPUTime: 0.20s
% Computational Cost: add. (115->38), mult. (344->68), div. (0->0), fcn. (260->6), ass. (0->36)
t26 = sin(pkin(7));
t27 = cos(pkin(7));
t39 = t26 ^ 2 + t27 ^ 2;
t48 = t39 * qJD(3);
t49 = 0.2e1 * t48;
t31 = cos(qJ(2));
t38 = pkin(1) * qJD(2);
t36 = t31 * t38;
t18 = qJD(3) + t36;
t47 = t39 * t18;
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t11 = t28 * t26 - t30 * t27;
t20 = -t27 * pkin(3) - pkin(2);
t42 = t31 * pkin(1);
t15 = t20 - t42;
t29 = sin(qJ(2));
t37 = t29 * t38;
t12 = t30 * t26 + t28 * t27;
t8 = t12 * qJD(4);
t46 = t11 * t37 + t15 * t8;
t7 = t11 * qJD(4);
t45 = t12 * t37 - t15 * t7;
t44 = t20 * t7;
t43 = t20 * t8;
t33 = t26 * t37;
t32 = t27 * t37;
t23 = t27 * pkin(6);
t19 = t29 * pkin(1) + qJ(3);
t17 = t27 * qJ(3) + t23;
t16 = (-pkin(6) - qJ(3)) * t26;
t10 = t27 * t19 + t23;
t9 = (-pkin(6) - t19) * t26;
t2 = -0.2e1 * t12 * t7;
t1 = 0.2e1 * t7 * t11 - 0.2e1 * t12 * t8;
t3 = [0, 0, 0, 0, -0.2e1 * t37, -0.2e1 * t36, -0.2e1 * t32, 0.2e1 * t33, 0.2e1 * t47, 0.2e1 * (-pkin(2) - t42) * t37 + 0.2e1 * t19 * t47, t2, t1, 0, 0, 0, 0.2e1 * t46, 0.2e1 * t45; 0, 0, 0, 0, -t37, -t36, -t32, t33, t48 + t47, -pkin(2) * t37 + qJ(3) * t47 + t19 * t48, t2, t1, 0, 0, 0, t43 + t46, -t44 + t45; 0, 0, 0, 0, 0, 0, 0, 0, t49, qJ(3) * t49, t2, t1, 0, 0, 0, 0.2e1 * t43, -0.2e1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, 0, 0, 0, 0, t8, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, -t12 * t18 + (-t10 * t30 - t28 * t9) * qJD(4), t11 * t18 + (t10 * t28 - t30 * t9) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, (-t16 * t28 - t17 * t30) * qJD(4) - t12 * qJD(3), (-t16 * t30 + t17 * t28) * qJD(4) + t11 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
