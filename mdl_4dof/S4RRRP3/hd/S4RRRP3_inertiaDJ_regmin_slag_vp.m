% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:15
% EndTime: 2019-12-31 17:14:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (122->41), mult. (343->74), div. (0->0), fcn. (192->4), ass. (0->36)
t39 = 2 * qJD(4);
t23 = sin(qJ(2));
t35 = pkin(1) * qJD(2);
t31 = t23 * t35;
t24 = cos(qJ(3));
t19 = qJD(3) * t24;
t22 = sin(qJ(3));
t34 = qJD(3) * t22;
t7 = pkin(3) * t34 - qJ(4) * t19 - t22 * qJD(4);
t1 = t7 + t31;
t38 = -t1 - t7;
t25 = cos(qJ(2));
t37 = t25 * pkin(1);
t17 = -pkin(2) - t37;
t36 = t17 * t19 + t22 * t31;
t33 = pkin(2) * t34;
t32 = pkin(2) * t19;
t30 = t25 * t35;
t29 = pkin(6) * t34;
t28 = pkin(6) * t19;
t27 = -t24 * pkin(3) - t22 * qJ(4);
t20 = t22 ^ 2;
t21 = t24 ^ 2;
t5 = (t20 + t21) * t30;
t26 = t17 * t34 - t24 * t31;
t13 = -pkin(2) + t27;
t6 = t27 * qJD(3) + t24 * qJD(4);
t16 = t23 * pkin(1) + pkin(6);
t15 = 0.2e1 * t22 * t19;
t10 = 0.2e1 * (-t20 + t21) * qJD(3);
t9 = t13 - t37;
t8 = t13 * t34;
t4 = t9 * t34;
t3 = t16 * t19 + t22 * t30;
t2 = t16 * t34 - t24 * t30;
t11 = [0, 0, 0, 0, -0.2e1 * t31, -0.2e1 * t30, t15, t10, 0, 0, 0, 0.2e1 * t26, 0.2e1 * t36, -0.2e1 * t1 * t24 + 0.2e1 * t4, 0.2e1 * t5, -0.2e1 * t1 * t22 - 0.2e1 * t9 * t19, 0.2e1 * t9 * t1 + 0.2e1 * t16 * t5; 0, 0, 0, 0, -t31, -t30, t15, t10, 0, 0, 0, t26 - t33, -t32 + t36, t38 * t24 + t4 + t8, t5, t38 * t22 + (-t13 - t9) * t19, pkin(6) * t5 + t1 * t13 + t9 * t7; 0, 0, 0, 0, 0, 0, t15, t10, 0, 0, 0, -0.2e1 * t33, -0.2e1 * t32, -0.2e1 * t7 * t24 + 0.2e1 * t8, 0, -0.2e1 * t13 * t19 - 0.2e1 * t7 * t22, 0.2e1 * t13 * t7; 0, 0, 0, 0, 0, 0, 0, 0, t19, -t34, 0, -t3, t2, -t3, t6, -t2, (-pkin(3) * t22 + qJ(4) * t24) * t30 + t6 * t16; 0, 0, 0, 0, 0, 0, 0, 0, t19, -t34, 0, -t28, t29, -t28, t6, -t29, t6 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, qJ(4) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
