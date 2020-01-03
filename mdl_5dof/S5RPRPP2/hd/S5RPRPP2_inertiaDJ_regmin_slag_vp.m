% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPP2
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
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:11:12
% EndTime: 2019-12-31 18:11:13
% DurationCPUTime: 0.19s
% Computational Cost: add. (134->45), mult. (281->78), div. (0->0), fcn. (177->4), ass. (0->31)
t17 = cos(qJ(3));
t32 = pkin(3) + pkin(4);
t33 = t32 * t17;
t19 = 2 * qJD(4);
t31 = t17 * pkin(3);
t9 = sin(pkin(7)) * pkin(1) + pkin(6);
t30 = qJ(5) - t9;
t29 = qJ(4) * t17;
t16 = sin(qJ(3));
t28 = t16 * qJ(4);
t12 = t16 * qJD(3);
t27 = t16 * qJD(4);
t13 = t17 * qJD(3);
t26 = t17 * qJD(4);
t25 = 0.2e1 * t13;
t24 = t9 * t12;
t23 = t9 * t13;
t10 = -cos(pkin(7)) * pkin(1) - pkin(2);
t8 = t30 * t17;
t22 = -pkin(3) * t12 + t27;
t21 = -t10 + t28;
t20 = t26 + (-t28 - t31) * qJD(3);
t14 = qJ(4) * t19;
t7 = t30 * t16;
t6 = -t21 - t31;
t5 = qJ(4) * t13 + t22;
t4 = t21 + t33;
t3 = -qJD(3) * t8 - t16 * qJD(5);
t2 = -t17 * qJD(5) + t30 * t12;
t1 = (-pkin(4) * t16 + t29) * qJD(3) + t22;
t11 = [0, 0, 0, 0, t16 * t25, 0.2e1 * (-t16 ^ 2 + t17 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t10 * t12, t10 * t25, 0.2e1 * t6 * t12 + 0.2e1 * t5 * t17, 0, -0.2e1 * t6 * t13 + 0.2e1 * t5 * t16, -0.2e1 * t6 * t5, 0.2e1 * t1 * t17 - 0.2e1 * t4 * t12, 0.2e1 * t1 * t16 + 0.2e1 * t4 * t13, -0.2e1 * t3 * t16 - 0.2e1 * t2 * t17 + 0.2e1 * (-t16 * t8 + t17 * t7) * qJD(3), 0.2e1 * t4 * t1 - 0.2e1 * t8 * t2 - 0.2e1 * t7 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2 * t16 - t3 * t17 + (-t16 * t7 - t17 * t8) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t13, -t12, 0, -t23, t24, -t23, t20, -t24, t20 * t9, -t3, t2, -t26 + (t28 + t33) * qJD(3), t2 * qJ(4) - t8 * qJD(4) - t3 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t13, -t12, 0, t13, t5, -t12, t13, 0, t27 + (-t16 * t32 + t29) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t14, 0, t19, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t23, 0, 0, -t13, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t13, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
