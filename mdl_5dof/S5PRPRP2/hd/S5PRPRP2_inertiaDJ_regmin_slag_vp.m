% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:04
% EndTime: 2019-12-05 15:31:05
% DurationCPUTime: 0.20s
% Computational Cost: add. (136->40), mult. (356->84), div. (0->0), fcn. (266->4), ass. (0->33)
t15 = sin(qJ(4));
t13 = sin(pkin(8));
t33 = qJ(5) * t13;
t14 = cos(pkin(8));
t9 = -t14 * pkin(3) - t13 * pkin(6) - pkin(2);
t21 = -t9 + t33;
t16 = cos(qJ(4));
t26 = t16 * t14 * qJ(3);
t38 = t21 * t15 - t26;
t28 = qJD(5) * t13;
t34 = qJ(3) * t15;
t29 = qJD(4) * t16;
t31 = qJD(3) * t16;
t35 = t14 * t31 + t9 * t29;
t1 = -t15 * t28 + (-t14 * t34 - t16 * t33) * qJD(4) + t35;
t32 = qJD(3) * t15;
t23 = t14 * t32;
t2 = t38 * qJD(4) - t16 * t28 - t23;
t3 = -t21 * t16 + (-pkin(4) - t34) * t14;
t37 = -t1 * t15 - t2 * t16 + (t15 * t3 + t16 * t38) * qJD(4);
t36 = 0.2e1 * qJD(4);
t30 = qJD(4) * t15;
t27 = qJ(3) * qJD(4);
t25 = t13 * t30;
t24 = t13 * t29;
t22 = t15 * t27;
t20 = t13 * t14 * t36;
t11 = t13 ^ 2;
t19 = 0.2e1 * (t14 ^ 2 + t11) * qJD(3);
t8 = (pkin(4) * t29 + qJD(3)) * t13;
t5 = -t23 + (-t15 * t9 - t26) * qJD(4);
t4 = t14 * t22 - t35;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t8 + (t1 * t16 - t15 * t2 + (t15 * t38 - t16 * t3) * qJD(4)) * t13; 0, 0, 0, 0, 0, 0, t19, qJ(3) * t19, -0.2e1 * t11 * t15 * t29, (t15 ^ 2 - t16 ^ 2) * t11 * t36, t15 * t20, t16 * t20, 0, -0.2e1 * t5 * t14 + 0.2e1 * (t16 * t27 + t32) * t11, -0.2e1 * t4 * t14 + 0.2e1 * (-t22 + t31) * t11, 0.2e1 * t37 * t13, -0.2e1 * t38 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * (pkin(4) * t15 + qJ(3)) * t8 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t30, t14 * t29, 0, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t25, 0, -pkin(4) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t24, 0, t5, t4, pkin(4) * t25, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t29, 0, -pkin(4) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
