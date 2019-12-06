% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:19
% EndTime: 2019-12-05 15:41:19
% DurationCPUTime: 0.17s
% Computational Cost: add. (74->42), mult. (197->73), div. (0->0), fcn. (130->4), ass. (0->28)
t27 = 2 * qJD(3);
t26 = 2 * qJD(5);
t12 = sin(qJ(4));
t10 = t12 ^ 2;
t14 = cos(qJ(4));
t11 = t14 ^ 2;
t25 = t10 + t11;
t24 = t12 * qJD(4);
t13 = sin(qJ(2));
t9 = t13 * qJD(2);
t23 = t14 * qJD(4);
t15 = cos(qJ(2));
t22 = t15 * qJD(2);
t21 = qJ(3) * qJD(4);
t16 = -pkin(2) - pkin(6);
t20 = t16 * t24;
t19 = t16 * t23;
t18 = t25 * t13;
t17 = -t12 * pkin(4) + t14 * qJ(5);
t2 = t17 * qJD(4) + t12 * qJD(5);
t8 = qJ(3) - t17;
t7 = -t13 * t24 + t14 * t22;
t6 = t12 * t22 + t13 * t23;
t5 = -t12 * t9 + t15 * t23;
t4 = t14 * t9 + t15 * t24;
t3 = qJD(2) * t18;
t1 = -t14 * qJD(5) + qJD(3) + (pkin(4) * t14 + qJ(5) * t12) * qJD(4);
t28 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (0.1e1 - t25) * t13 * t22; 0, 0, -t9, -t22, t9, t22, t13 * qJD(3) + (-pkin(2) * t13 + qJ(3) * t15) * qJD(2), 0, 0, 0, 0, 0, t6, t7, t6, -t3, -t7, t13 * t1 + (t15 * t8 + t16 * t18) * qJD(2); 0, 0, 0, 0, 0, t27, qJ(3) * t27, -0.2e1 * t12 * t23, 0.2e1 * (t10 - t11) * qJD(4), 0, 0, 0, 0.2e1 * qJD(3) * t12 + 0.2e1 * t14 * t21, 0.2e1 * qJD(3) * t14 - 0.2e1 * t12 * t21, 0.2e1 * t1 * t12 + 0.2e1 * t8 * t23, 0, -0.2e1 * t1 * t14 + 0.2e1 * t8 * t24, 0.2e1 * t8 * t1; 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, t4, 0, -t5, (-qJ(5) * qJD(4) * t15 + pkin(4) * t9) * t14 + (qJ(5) * t9 + (qJD(4) * pkin(4) - qJD(5)) * t15) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, 0, -t20, -t19, -t20, -t2, t19, t2 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, -t24, 0, t23, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, qJ(5) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t28;
