% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPPP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:26
% EndTime: 2019-03-08 18:26:26
% DurationCPUTime: 0.14s
% Computational Cost: add. (46->29), mult. (180->80), div. (0->0), fcn. (139->4), ass. (0->28)
t11 = sin(pkin(6));
t13 = cos(pkin(6));
t12 = sin(pkin(4));
t25 = qJ(2) * t12;
t14 = cos(pkin(4));
t28 = pkin(1) * t14;
t29 = t11 * t28 + t13 * t25;
t23 = qJD(2) * t12;
t3 = t14 * qJD(3) + t13 * t23;
t27 = t3 * t14;
t26 = t12 * t13;
t10 = t12 ^ 2;
t24 = qJD(2) * t10;
t22 = qJD(3) * t10;
t21 = t11 * qJD(3);
t20 = t14 * t23;
t19 = t12 * t21;
t18 = -pkin(1) * t13 - pkin(2);
t17 = -qJ(3) * t11 - pkin(1);
t16 = -0.2e1 * t20;
t5 = t11 * t23;
t2 = -t14 * qJD(4) + t5;
t15 = t11 * t2 + t13 * t3;
t9 = t11 ^ 2;
t6 = t11 * t25;
t4 = t9 * t24;
t1 = (-qJD(4) * t13 - t21) * t12;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t16, t13 * t16, 0.2e1 * t13 ^ 2 * t24 + 0.2e1 * t4, 0.2e1 * (t29 * t13 - (t13 * t28 - t6) * t11) * t23, 0, 0, 0, 0, 0, 0, 0.2e1 * t3 * t26 + 0.2e1 * t4, 0.2e1 * (-t13 * t22 + t20) * t11, 0.2e1 * t9 * t22 + 0.2e1 * t27, -0.2e1 * (-t14 * qJ(3) - t29) * t3 + 0.2e1 * ((t18 * t14 + t6) * t11 * qJD(2) - (-pkin(2) * t13 + t17) * t19) * t12, 0, 0, 0, 0, 0, 0, 0.2e1 * t15 * t12, -0.2e1 * t1 * t12 * t11 + 0.2e1 * t27, -0.2e1 * t1 * t26 - 0.2e1 * t2 * t14, 0.2e1 * t6 * t2 + 0.2e1 * t29 * t3 + 0.2e1 * ((-qJ(4) + t18) * t2 + qJ(3) * t3) * t14 + 0.2e1 * (((-pkin(2) - qJ(4)) * t13 + t17) * t1 + t15 * pkin(3)) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
