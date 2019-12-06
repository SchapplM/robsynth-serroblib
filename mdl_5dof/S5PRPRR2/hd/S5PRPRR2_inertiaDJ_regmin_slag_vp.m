% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x15]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:13
% EndTime: 2019-12-05 15:45:14
% DurationCPUTime: 0.18s
% Computational Cost: add. (158->38), mult. (421->69), div. (0->0), fcn. (400->8), ass. (0->37)
t28 = cos(pkin(9));
t25 = t28 * pkin(2) + pkin(3);
t30 = sin(qJ(4));
t33 = cos(qJ(4));
t27 = sin(pkin(9));
t43 = pkin(2) * t27;
t35 = t30 * t25 + t33 * t43;
t11 = t35 * qJD(4);
t39 = t30 * t43;
t15 = -t33 * t25 - pkin(4) + t39;
t32 = cos(qJ(5));
t26 = qJD(5) * t32;
t29 = sin(qJ(5));
t44 = t11 * t29 + t15 * t26;
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t21 = t27 * t34 + t28 * t31;
t42 = t30 * t21;
t41 = qJD(4) * t33;
t40 = qJD(5) * t29;
t38 = pkin(4) * t40;
t37 = pkin(4) * t26;
t36 = -t11 * t32 + t15 * t40;
t20 = -t27 * t31 + t28 * t34;
t6 = t30 * t20 + t33 * t21;
t24 = 0.2e1 * t29 * t26;
t22 = 0.2e1 * (-t29 ^ 2 + t32 ^ 2) * qJD(5);
t19 = t20 * qJD(2);
t18 = t21 * qJD(2);
t16 = pkin(7) + t35;
t10 = qJD(4) * t39 - t25 * t41;
t5 = -t33 * t20 + t42;
t4 = t6 * qJD(4) + t33 * t18 + t30 * t19;
t3 = qJD(4) * t42 + t30 * t18 - t33 * t19 - t20 * t41;
t2 = -t4 * t32 + t5 * t40;
t1 = t5 * t26 + t4 * t29;
t7 = [0, 0, 0, 0, -0.2e1 * t20 * t18 + 0.2e1 * t21 * t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -qJD(2) * t31, -qJD(2) * t34, (-t18 * t28 + t19 * t27) * pkin(2), 0, -t4, t3, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, -0.2e1 * t11, 0.2e1 * t10, t24, t22, 0, 0, 0, 0.2e1 * t36, 0.2e1 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -t4, t3, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, -t11, t10, t24, t22, 0, 0, 0, t36 - t38, -t37 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t24, t22, 0, 0, 0, -0.2e1 * t38, -0.2e1 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t26 + t29 * t3, t32 * t3 + t6 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t40, 0, t29 * t10 - t16 * t26, t32 * t10 + t16 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t40, 0, -pkin(7) * t26, pkin(7) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
