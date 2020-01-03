% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPR3
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
% MMD_reg [((4+1)*4/2)x14]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:34
% DurationCPUTime: 0.12s
% Computational Cost: add. (58->25), mult. (185->48), div. (0->0), fcn. (120->6), ass. (0->29)
t20 = cos(qJ(4));
t15 = t20 * qJD(4);
t18 = sin(qJ(4));
t21 = cos(qJ(2));
t14 = t21 * pkin(1) + pkin(2);
t17 = cos(pkin(7));
t16 = sin(pkin(7));
t19 = sin(qJ(2));
t31 = t16 * t19;
t22 = -pkin(1) * t31 + t17 * t14;
t4 = -pkin(3) - t22;
t29 = pkin(1) * qJD(2);
t30 = t17 * t19;
t6 = (t16 * t21 + t30) * t29;
t33 = t4 * t15 + t6 * t18;
t32 = pkin(1) * t30 + t16 * t14;
t28 = t18 * qJD(4);
t27 = t19 * t29;
t26 = t21 * t29;
t13 = -t17 * pkin(2) - pkin(3);
t25 = t13 * t28;
t24 = t13 * t15;
t23 = -t6 * t20 + t4 * t28;
t12 = t16 * pkin(2) + pkin(6);
t10 = 0.2e1 * t18 * t15;
t8 = 0.2e1 * (-t18 ^ 2 + t20 ^ 2) * qJD(4);
t7 = (t17 * t21 - t31) * t29;
t5 = pkin(6) + t32;
t1 = [0, 0, 0, 0, -0.2e1 * t27, -0.2e1 * t26, -0.2e1 * t22 * t6 + 0.2e1 * t32 * t7, t10, t8, 0, 0, 0, 0.2e1 * t23, 0.2e1 * t33; 0, 0, 0, 0, -t27, -t26, (t16 * t7 - t17 * t6) * pkin(2), t10, t8, 0, 0, 0, t23 + t25, t24 + t33; 0, 0, 0, 0, 0, 0, 0, t10, t8, 0, 0, 0, 0.2e1 * t25, 0.2e1 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t28, 0, -t5 * t15 - t18 * t7, -t20 * t7 + t5 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t28, 0, -t12 * t15, t12 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
