% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4PRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:43
% EndTime: 2019-12-31 16:30:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (48->27), mult. (166->55), div. (0->0), fcn. (110->4), ass. (0->26)
t25 = 2 * qJD(4);
t10 = sin(qJ(3));
t8 = t10 ^ 2;
t12 = cos(qJ(3));
t9 = t12 ^ 2;
t24 = t8 + t9;
t23 = t10 * qJD(3);
t11 = sin(qJ(2));
t22 = t11 * qJD(2);
t7 = t12 * qJD(3);
t13 = cos(qJ(2));
t21 = t13 * qJD(2);
t20 = -0.2e1 * pkin(2) * qJD(3);
t19 = pkin(5) * t23;
t18 = pkin(5) * t7;
t17 = t24 * t13;
t16 = -t12 * pkin(3) - t10 * qJ(4);
t15 = pkin(3) * t10 - qJ(4) * t12;
t14 = t16 * qJD(3) + t12 * qJD(4);
t6 = -pkin(2) + t16;
t5 = -t12 * t22 - t13 * t23;
t4 = t10 * t22 - t13 * t7;
t3 = t10 * t21 + t11 * t7;
t2 = t11 * t23 - t12 * t21;
t1 = t15 * qJD(3) - t10 * qJD(4);
t26 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t24) * t11 * t21; 0, 0, -t22, -t21, 0, 0, 0, 0, 0, t5, t4, t5, qJD(2) * t17, -t4, -t13 * t1 + (pkin(5) * t17 + t11 * t6) * qJD(2); 0, 0, 0, 0, 0.2e1 * t10 * t7, 0.2e1 * (-t8 + t9) * qJD(3), 0, 0, 0, t10 * t20, t12 * t20, -0.2e1 * t1 * t12 + 0.2e1 * t6 * t23, 0, -0.2e1 * t1 * t10 - 0.2e1 * t6 * t7, 0.2e1 * t6 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, t2, -t3, 0, -t2, t14 * t11 - t15 * t21; 0, 0, 0, 0, 0, 0, t7, -t23, 0, -t18, t19, -t18, t14, -t19, t14 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, qJ(4) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t26;
