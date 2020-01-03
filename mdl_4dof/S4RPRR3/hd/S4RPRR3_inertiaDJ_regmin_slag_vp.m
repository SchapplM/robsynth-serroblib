% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x18]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:18
% EndTime: 2019-12-31 16:49:19
% DurationCPUTime: 0.14s
% Computational Cost: add. (106->31), mult. (261->59), div. (0->0), fcn. (214->6), ass. (0->29)
t31 = qJD(3) + qJD(4);
t14 = sin(pkin(7)) * pkin(1) + pkin(5);
t30 = pkin(6) + t14;
t17 = sin(qJ(4));
t18 = sin(qJ(3));
t29 = t17 * t18;
t28 = qJD(3) * t18;
t20 = cos(qJ(3));
t27 = qJD(3) * t20;
t19 = cos(qJ(4));
t26 = qJD(4) * t19;
t25 = 0.2e1 * t27;
t24 = pkin(3) * t28;
t23 = qJD(4) * t17 * pkin(3);
t22 = pkin(3) * t26;
t15 = -cos(pkin(7)) * pkin(1) - pkin(2);
t21 = qJD(3) * t30;
t10 = t17 * t20 + t19 * t18;
t11 = -t20 * pkin(3) + t15;
t9 = -t19 * t20 + t29;
t8 = t30 * t20;
t7 = t30 * t18;
t6 = t20 * t21;
t5 = t18 * t21;
t4 = t31 * t10;
t3 = -t19 * t27 - t20 * t26 + t31 * t29;
t2 = t17 * t5 - t19 * t6 + (t17 * t7 - t19 * t8) * qJD(4);
t1 = t17 * t6 + t19 * t5 + (t17 * t8 + t19 * t7) * qJD(4);
t12 = [0, 0, 0, 0, t18 * t25, 0.2e1 * (-t18 ^ 2 + t20 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * t15 * t28, t15 * t25, -0.2e1 * t10 * t3, -0.2e1 * t10 * t4 + 0.2e1 * t3 * t9, 0, 0, 0, 0.2e1 * t11 * t4 + 0.2e1 * t9 * t24, 0.2e1 * t10 * t24 - 0.2e1 * t11 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t27, -t28, 0, -t14 * t27, t14 * t28, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t27, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t23, -0.2e1 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t12;
