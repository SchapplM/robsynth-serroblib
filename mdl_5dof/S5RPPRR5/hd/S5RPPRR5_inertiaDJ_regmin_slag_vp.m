% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:38
% EndTime: 2019-12-31 17:56:38
% DurationCPUTime: 0.14s
% Computational Cost: add. (103->30), mult. (201->56), div. (0->0), fcn. (139->6), ass. (0->28)
t29 = 2 * qJD(3);
t16 = sin(qJ(5));
t17 = sin(qJ(4));
t11 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(3);
t12 = sin(pkin(8)) * pkin(1) + qJ(3);
t19 = cos(qJ(4));
t20 = t17 * t11 + t19 * t12;
t2 = t17 * qJD(3) + t20 * qJD(4);
t28 = t2 * t16;
t18 = cos(qJ(5));
t27 = t2 * t18;
t26 = qJD(4) * t17;
t25 = qJD(4) * t19;
t14 = qJD(5) * t16;
t15 = qJD(5) * t18;
t24 = qJD(5) * t19;
t23 = -0.2e1 * pkin(4) * qJD(5);
t22 = t16 * t15;
t3 = -t19 * t11 + t17 * t12 + pkin(4);
t21 = qJD(5) * (pkin(4) + t3);
t10 = 0.2e1 * t22;
t9 = (-t16 ^ 2 + t18 ^ 2) * qJD(5);
t7 = 0.2e1 * t9;
t6 = t16 * t24 + t18 * t26;
t5 = t16 * t26 - t18 * t24;
t4 = -pkin(7) + t20;
t1 = -t19 * qJD(3) - t11 * t25 + t12 * t26;
t8 = [0, 0, 0, 0, 0, t29, t12 * t29, 0, 0.2e1 * t2, -0.2e1 * t1, t10, t7, 0, 0, 0, -0.2e1 * t3 * t14 + 0.2e1 * t27, -0.2e1 * t3 * t15 - 0.2e1 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t26, t25, 0, 0, 0, 0, 0, t6, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -0.2e1 * t22, -0.2e1 * t9, 0, 0, 0, t16 * t21 - t27, t18 * t21 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t26, -t25, 0, 0, 0, 0, 0, -t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t7, 0, 0, 0, t16 * t23, t18 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t14, 0, t16 * t1 - t4 * t15, t18 * t1 + t4 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * t15 - t16 * t25, t17 * t14 - t18 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, -pkin(7) * t15, pkin(7) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
