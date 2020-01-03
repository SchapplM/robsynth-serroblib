% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPP5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:42
% EndTime: 2019-12-31 18:16:42
% DurationCPUTime: 0.20s
% Computational Cost: add. (137->45), mult. (254->75), div. (0->0), fcn. (145->2), ass. (0->31)
t16 = sin(qJ(3));
t18 = -pkin(3) - pkin(4);
t17 = cos(qJ(3));
t27 = t17 * qJ(4);
t30 = t16 * t18 + t27;
t29 = 2 * qJD(2);
t20 = 2 * qJD(4);
t28 = qJ(4) * t16;
t19 = -pkin(1) - pkin(6);
t26 = qJ(5) + t19;
t13 = t16 * qJD(3);
t25 = t16 * qJD(4);
t14 = t17 * qJD(3);
t24 = qJ(2) * qJD(3);
t23 = t19 * t13;
t9 = t26 * t16;
t22 = -t17 * qJD(4) + qJD(2);
t21 = -t16 * pkin(3) + t27;
t7 = t21 * qJD(3) + t25;
t15 = qJ(4) * t20;
t12 = t19 * t14;
t11 = qJ(2) - t21;
t10 = t26 * t17;
t8 = -qJ(2) + t30;
t6 = (pkin(3) * t17 + t28) * qJD(3) + t22;
t5 = t30 * qJD(3) + t25;
t4 = qJ(5) * t14 + t16 * qJD(5) + t12;
t3 = qJD(3) * t9 - t17 * qJD(5);
t2 = (t18 * t17 - t28) * qJD(3) - t22;
t1 = t4 * t16 - t3 * t17 + (-t10 * t16 + t17 * t9) * qJD(3);
t31 = [0, 0, 0, 0, t29, qJ(2) * t29, -0.2e1 * t16 * t14, 0.2e1 * (t16 ^ 2 - t17 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t16 + 0.2e1 * t17 * t24, 0.2e1 * qJD(2) * t17 - 0.2e1 * t16 * t24, 0.2e1 * t11 * t14 + 0.2e1 * t6 * t16, 0, 0.2e1 * t11 * t13 - 0.2e1 * t6 * t17, 0.2e1 * t11 * t6, -0.2e1 * t8 * t14 - 0.2e1 * t2 * t16, -0.2e1 * t8 * t13 + 0.2e1 * t2 * t17, 0.2e1 * t1, -0.2e1 * t10 * t3 + 0.2e1 * t8 * t2 + 0.2e1 * t9 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, 0, -t23, -t12, -t23, -t7, t12, t7 * t19, -t3, t4, t5, t4 * qJ(4) + t9 * qJD(4) + t3 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, -t13, 0, t14, t7, -t13, t14, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t15, 0, t20, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, t23, 0, 0, t13, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t13, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t31;
