% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:45
% EndTime: 2019-12-31 17:39:45
% DurationCPUTime: 0.14s
% Computational Cost: add. (70->28), mult. (168->54), div. (0->0), fcn. (106->4), ass. (0->27)
t28 = 2 * qJD(3);
t14 = sin(qJ(5));
t15 = sin(qJ(4));
t17 = cos(qJ(4));
t18 = -pkin(2) - pkin(3);
t19 = t17 * qJ(3) + t15 * t18;
t2 = t15 * qJD(3) + t19 * qJD(4);
t27 = t2 * t14;
t16 = cos(qJ(5));
t26 = t2 * t16;
t25 = qJD(4) * t15;
t24 = qJD(4) * t17;
t12 = qJD(5) * t14;
t13 = qJD(5) * t16;
t23 = qJD(5) * t17;
t22 = -0.2e1 * pkin(4) * qJD(5);
t21 = t14 * t13;
t6 = t15 * qJ(3) - t17 * t18 + pkin(4);
t20 = qJD(5) * (pkin(4) + t6);
t9 = 0.2e1 * t21;
t8 = (-t14 ^ 2 + t16 ^ 2) * qJD(5);
t7 = -pkin(7) + t19;
t5 = 0.2e1 * t8;
t4 = t14 * t23 + t16 * t25;
t3 = t14 * t25 - t16 * t23;
t1 = qJ(3) * t25 - t17 * qJD(3) - t18 * t24;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t28, qJ(3) * t28, 0, 0.2e1 * t2, -0.2e1 * t1, t9, t5, 0, 0, 0, -0.2e1 * t6 * t12 + 0.2e1 * t26, -0.2e1 * t6 * t13 - 0.2e1 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t25, t24, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -0.2e1 * t21, -0.2e1 * t8, 0, 0, 0, t14 * t20 - t26, t16 * t20 + t27; 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t24, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t5, 0, 0, 0, t14 * t22, t16 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t12, 0, t14 * t1 - t7 * t13, t16 * t1 + t7 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t13 - t14 * t24, t15 * t12 - t16 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, 0, -pkin(7) * t13, pkin(7) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
