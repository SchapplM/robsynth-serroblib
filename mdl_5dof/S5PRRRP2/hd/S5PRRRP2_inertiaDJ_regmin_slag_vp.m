% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:57
% EndTime: 2019-12-05 16:41:59
% DurationCPUTime: 0.20s
% Computational Cost: add. (127->44), mult. (353->74), div. (0->0), fcn. (200->4), ass. (0->36)
t39 = 2 * qJD(5);
t23 = sin(qJ(3));
t35 = pkin(2) * qJD(3);
t31 = t23 * t35;
t24 = cos(qJ(4));
t19 = qJD(4) * t24;
t22 = sin(qJ(4));
t34 = qJD(4) * t22;
t7 = -pkin(4) * t34 + qJ(5) * t19 + t22 * qJD(5);
t1 = -t7 + t31;
t38 = -t1 + t7;
t25 = cos(qJ(3));
t37 = t25 * pkin(2);
t17 = -pkin(3) - t37;
t36 = t17 * t19 + t22 * t31;
t33 = pkin(3) * t34;
t32 = pkin(3) * t19;
t30 = t25 * t35;
t29 = pkin(7) * t34;
t28 = pkin(7) * t19;
t27 = -t24 * pkin(4) - t22 * qJ(5);
t20 = t22 ^ 2;
t21 = t24 ^ 2;
t5 = (t20 + t21) * t30;
t26 = t17 * t34 - t24 * t31;
t13 = -pkin(3) + t27;
t6 = t27 * qJD(4) + t24 * qJD(5);
t16 = t23 * pkin(2) + pkin(7);
t15 = 0.2e1 * t22 * t19;
t10 = 0.2e1 * (-t20 + t21) * qJD(4);
t9 = t13 - t37;
t8 = t13 * t34;
t4 = t9 * t34;
t3 = t16 * t19 + t22 * t30;
t2 = t16 * t34 - t24 * t30;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -0.2e1 * t31, -0.2e1 * t30, t15, t10, 0, 0, 0, 0.2e1 * t26, 0.2e1 * t36, -0.2e1 * t1 * t24 + 0.2e1 * t4, 0.2e1 * t5, -0.2e1 * t1 * t22 - 0.2e1 * t9 * t19, 0.2e1 * t9 * t1 + 0.2e1 * t16 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t31, -t30, t15, t10, 0, 0, 0, t26 - t33, -t32 + t36, t38 * t24 + t4 + t8, t5, t38 * t22 + (-t13 - t9) * t19, pkin(7) * t5 + t1 * t13 - t9 * t7; 0, 0, 0, 0, 0, 0, 0, t15, t10, 0, 0, 0, -0.2e1 * t33, -0.2e1 * t32, 0.2e1 * t7 * t24 + 0.2e1 * t8, 0, -0.2e1 * t13 * t19 + 0.2e1 * t7 * t22, -0.2e1 * t13 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t19, -t34, 0, t19, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t34, 0, -t3, t2, -t3, t6, -t2, (-pkin(4) * t22 + qJ(5) * t24) * t30 + t6 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t34, 0, -t28, t29, -t28, t6, -t29, t6 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, qJ(5) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
