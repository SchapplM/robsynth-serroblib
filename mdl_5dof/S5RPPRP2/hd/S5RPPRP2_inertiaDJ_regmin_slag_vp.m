% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:24
% EndTime: 2019-12-31 17:49:24
% DurationCPUTime: 0.18s
% Computational Cost: add. (243->47), mult. (556->84), div. (0->0), fcn. (500->6), ass. (0->29)
t37 = 2 * qJD(5);
t18 = sin(pkin(7)) * pkin(1) + qJ(3);
t36 = pkin(6) + t18;
t35 = cos(qJ(4));
t23 = cos(pkin(8));
t29 = qJD(4) * t35;
t22 = sin(pkin(8));
t25 = sin(qJ(4));
t33 = t25 * t22;
t10 = qJD(4) * t33 - t23 * t29;
t30 = t35 * t22;
t32 = t25 * t23;
t14 = t30 + t32;
t34 = t14 * t10;
t31 = t25 * t36;
t28 = t35 * qJD(3);
t27 = 0.2e1 * (t22 ^ 2 + t23 ^ 2) * qJD(3);
t15 = -cos(pkin(7)) * pkin(1) - pkin(2) - t23 * pkin(3);
t26 = t36 * t30;
t11 = t14 * qJD(4);
t3 = t11 * pkin(4) + t10 * qJ(5) - t14 * qJD(5);
t13 = -t35 * t23 + t33;
t12 = t36 * t23;
t6 = t35 * t12 - t22 * t31;
t5 = t25 * t12 + t26;
t4 = t13 * pkin(4) - t14 * qJ(5) + t15;
t2 = t12 * t29 + qJD(3) * t32 + (-qJD(4) * t31 + t28) * t22;
t1 = -t23 * t28 + qJD(4) * t26 + (qJD(3) * t22 + qJD(4) * t12) * t25;
t7 = [0, 0, 0, 0, 0, 0, t27, t18 * t27, -0.2e1 * t34, 0.2e1 * t13 * t10 - 0.2e1 * t11 * t14, 0, 0, 0, 0.2e1 * t15 * t11, -0.2e1 * t15 * t10, 0.2e1 * t4 * t11 + 0.2e1 * t3 * t13, 0.2e1 * t1 * t13 - 0.2e1 * t5 * t10 - 0.2e1 * t6 * t11 + 0.2e1 * t2 * t14, 0.2e1 * t4 * t10 - 0.2e1 * t3 * t14, -0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * t4 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t14 - t6 * t10 + t5 * t11 + t2 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t13 * t11 - 0.2e1 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, t11, 0, t10, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, -t2, t1, -t2, pkin(4) * t10 - t11 * qJ(5) - t13 * qJD(5), -t1, -t2 * pkin(4) - t1 * qJ(5) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, -t11, 0, -t10, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, qJ(5) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
