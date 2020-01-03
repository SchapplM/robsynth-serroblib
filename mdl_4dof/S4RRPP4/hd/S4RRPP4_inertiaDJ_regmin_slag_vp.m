% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPP4
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:18
% EndTime: 2019-12-31 16:59:19
% DurationCPUTime: 0.14s
% Computational Cost: add. (82->33), mult. (209->69), div. (0->0), fcn. (115->2), ass. (0->27)
t12 = cos(qJ(2));
t11 = sin(qJ(2));
t23 = t11 * qJ(3);
t26 = pkin(2) + pkin(3);
t27 = t26 * t12 + t23;
t14 = 2 * qJD(3);
t25 = pkin(5) - qJ(4);
t24 = qJ(3) * t12;
t22 = qJD(2) * t11;
t9 = qJD(2) * t12;
t21 = t11 * qJD(3);
t20 = t12 * qJD(3);
t19 = -0.2e1 * pkin(1) * qJD(2);
t18 = pkin(5) * t22;
t17 = pkin(5) * t9;
t8 = t25 * t12;
t16 = -t12 * pkin(2) - t23;
t15 = t16 * qJD(2) + t20;
t10 = qJ(3) * t14;
t7 = t25 * t11;
t6 = -pkin(1) + t16;
t5 = pkin(1) + t27;
t4 = -t21 + (pkin(2) * t11 - t24) * qJD(2);
t3 = qJD(2) * t8 - t11 * qJD(4);
t2 = -t12 * qJD(4) - t25 * t22;
t1 = t21 + (-t26 * t11 + t24) * qJD(2);
t13 = [0, 0, 0, 0.2e1 * t11 * t9, 0.2e1 * (-t11 ^ 2 + t12 ^ 2) * qJD(2), 0, 0, 0, t11 * t19, t12 * t19, -0.2e1 * t4 * t12 + 0.2e1 * t6 * t22, 0, -0.2e1 * t4 * t11 - 0.2e1 * t6 * t9, 0.2e1 * t6 * t4, 0.2e1 * t1 * t12 - 0.2e1 * t5 * t22, 0.2e1 * t1 * t11 + 0.2e1 * t5 * t9, -0.2e1 * t3 * t11 - 0.2e1 * t2 * t12 + 0.2e1 * (t11 * t8 - t12 * t7) * qJD(2), 0.2e1 * t5 * t1 + 0.2e1 * t8 * t2 + 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, t9, -t22, 0, -t17, t18, -t17, t15, -t18, t15 * pkin(5), -t3, t2, t27 * qJD(2) - t20, t2 * qJ(3) + t8 * qJD(3) - t26 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t10, 0, t14, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, t17, 0, 0, -t9, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t9, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
