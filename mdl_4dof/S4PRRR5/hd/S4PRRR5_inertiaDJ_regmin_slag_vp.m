% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x14]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:46
% EndTime: 2019-12-31 16:33:46
% DurationCPUTime: 0.13s
% Computational Cost: add. (66->27), mult. (207->48), div. (0->0), fcn. (168->6), ass. (0->30)
t34 = qJD(2) + qJD(3);
t18 = sin(qJ(3));
t33 = t18 * pkin(2);
t21 = cos(qJ(3));
t15 = -t21 * pkin(2) - pkin(3);
t20 = cos(qJ(4));
t16 = t20 * qJD(4);
t17 = sin(qJ(4));
t25 = qJD(3) * t33;
t32 = t15 * t16 + t17 * t25;
t19 = sin(qJ(2));
t31 = t18 * t19;
t30 = qJD(3) * t21;
t29 = t17 * qJD(4);
t22 = cos(qJ(2));
t28 = t22 * qJD(2);
t27 = pkin(3) * t29;
t26 = pkin(3) * t16;
t24 = pkin(2) * t30;
t7 = t18 * t22 + t21 * t19;
t23 = t15 * t29 - t20 * t25;
t14 = pkin(6) + t33;
t11 = 0.2e1 * t17 * t16;
t6 = -t21 * t22 + t31;
t5 = 0.2e1 * (-t17 ^ 2 + t20 ^ 2) * qJD(4);
t4 = t34 * t7;
t3 = -t21 * t28 - t22 * t30 + t34 * t31;
t2 = -t4 * t20 + t6 * t29;
t1 = t6 * t16 + t4 * t17;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t19 * qJD(2), -t28, 0, -t4, t3, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, -0.2e1 * t25, -0.2e1 * t24, t11, t5, 0, 0, 0, 0.2e1 * t23, 0.2e1 * t32; 0, 0, 0, 0, 0, -t4, t3, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, -t25, -t24, t11, t5, 0, 0, 0, t23 - t27, -t26 + t32; 0, 0, 0, 0, 0, 0, 0, t11, t5, 0, 0, 0, -0.2e1 * t27, -0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7 * t16 + t17 * t3, t20 * t3 + t7 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t29, 0, -t14 * t16 - t17 * t24, t14 * t29 - t20 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t29, 0, -pkin(6) * t16, pkin(6) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
