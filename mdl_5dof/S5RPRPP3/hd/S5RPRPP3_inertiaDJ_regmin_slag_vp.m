% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:54
% EndTime: 2019-12-31 18:12:55
% DurationCPUTime: 0.27s
% Computational Cost: add. (338->60), mult. (806->111), div. (0->0), fcn. (708->4), ass. (0->36)
t27 = sin(pkin(7));
t42 = pkin(6) + qJ(2);
t19 = t42 * t27;
t28 = cos(pkin(7));
t20 = t42 * t28;
t30 = sin(qJ(3));
t45 = cos(qJ(3));
t38 = qJD(3) * t45;
t39 = t45 * t28;
t6 = (qJD(2) * t27 + qJD(3) * t20) * t30 - qJD(2) * t39 + t19 * t38;
t44 = t30 * t27;
t13 = qJD(3) * t44 - t28 * t38;
t46 = -0.2e1 * t13;
t43 = pkin(3) + qJ(5);
t40 = (qJ(4) * qJD(4));
t24 = -t28 * pkin(2) - pkin(1);
t37 = 0.2e1 * (t27 ^ 2 + t28 ^ 2) * qJD(2);
t11 = t45 * t19 + t30 * t20;
t18 = t45 * t27 + t30 * t28;
t36 = t13 * qJ(4) - t18 * qJD(4);
t14 = t18 * qJD(3);
t17 = -t39 + t44;
t35 = -qJ(4) * t14 - qJD(4) * t17;
t33 = -t18 * qJ(4) + t24;
t32 = -t30 * t19 + t45 * t20;
t7 = t18 * qJD(2) + t32 * qJD(3);
t31 = 0.2e1 * qJD(4);
t10 = t17 * pkin(3) + t33;
t9 = -t17 * pkin(4) + t32;
t8 = t18 * pkin(4) + t11;
t5 = t43 * t17 + t33;
t4 = t14 * pkin(3) + t36;
t3 = -t13 * pkin(4) + t7;
t2 = -t14 * pkin(4) - t6;
t1 = t17 * qJD(5) + t43 * t14 + t36;
t12 = [0, 0, 0, 0, 0, t37, qJ(2) * t37, t18 * t46, 0.2e1 * t13 * t17 - 0.2e1 * t18 * t14, 0, 0, 0, 0.2e1 * t24 * t14, t24 * t46, -0.2e1 * t11 * t13 - 0.2e1 * t14 * t32 + 0.2e1 * t6 * t17 + 0.2e1 * t7 * t18, -0.2e1 * t10 * t14 - 0.2e1 * t4 * t17, 0.2e1 * t10 * t13 - 0.2e1 * t4 * t18, 0.2e1 * t10 * t4 + 0.2e1 * t11 * t7 - 0.2e1 * t32 * t6, -0.2e1 * t8 * t13 - 0.2e1 * t9 * t14 - 0.2e1 * t2 * t17 + 0.2e1 * t3 * t18, -0.2e1 * t1 * t18 + 0.2e1 * t5 * t13, 0.2e1 * t1 * t17 + 0.2e1 * t5 * t14, 0.2e1 * t5 * t1 + 0.2e1 * t9 * t2 + 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, -t14, t13, t4, 0, t13, t14, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, -t7, t6, pkin(3) * t13 + t35, t7, -t6, -t7 * pkin(3) - t6 * qJ(4) + qJD(4) * t32, -qJD(5) * t18 + t13 * t43 + t35, t2, -t3, t2 * qJ(4) + t9 * qJD(4) - t8 * qJD(5) - t3 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 2 * t40, 0, t31, 0.2e1 * qJD(5), 0.2e1 * qJD(5) * t43 + (2 * t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, t7, -t13, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
