% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:55
% EndTime: 2019-12-31 17:36:56
% DurationCPUTime: 0.13s
% Computational Cost: add. (45->24), mult. (143->51), div. (0->0), fcn. (110->4), ass. (0->21)
t14 = sin(pkin(8));
t12 = t14 ^ 2;
t15 = cos(pkin(8));
t25 = (t15 ^ 2 + t12) * qJD(3);
t24 = -pkin(6) + qJ(3);
t23 = qJ(3) * t25;
t16 = sin(qJ(5));
t22 = qJD(5) * t16;
t17 = cos(qJ(5));
t21 = qJD(5) * t17;
t20 = t14 * qJD(4);
t18 = t14 * qJ(4) + pkin(2);
t6 = t14 * t17 - t15 * t16;
t5 = t14 * t16 + t15 * t17;
t8 = t24 * t15;
t7 = t24 * t14;
t4 = 0.2e1 * t25;
t3 = (pkin(3) + pkin(4)) * t15 + t18;
t2 = t14 * t21 - t15 * t22;
t1 = t5 * qJD(5);
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t4, 0.2e1 * t23, 0.2e1 * t15 * t20, t4, 0.2e1 * t12 * qJD(4), -0.2e1 * (-t15 * pkin(3) - t18) * t20 + 0.2e1 * t23, -0.2e1 * t6 * t1, 0.2e1 * t1 * t5 - 0.2e1 * t6 * t2, 0, 0, 0, 0.2e1 * t3 * t2 + 0.2e1 * t5 * t20, -0.2e1 * t3 * t1 + 0.2e1 * t6 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * qJD(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, 0, (-t16 * t7 - t17 * t8) * qJD(5) + t6 * qJD(3), (t16 * t8 - t17 * t7) * qJD(5) - t5 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
