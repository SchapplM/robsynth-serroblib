% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x12]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:34
% EndTime: 2019-07-18 18:16:35
% DurationCPUTime: 0.12s
% Computational Cost: add. (77->34), mult. (167->44), div. (0->0), fcn. (88->4), ass. (0->25)
t22 = 2 * qJD(3);
t19 = cos(qJ(4));
t16 = qJD(4) * t19;
t21 = -pkin(2) - pkin(3);
t29 = t19 * qJD(3) + t21 * t16;
t28 = pkin(1) * qJD(2);
t18 = sin(qJ(2));
t12 = t18 * pkin(1) + qJ(3);
t27 = qJ(3) + t12;
t17 = sin(qJ(4));
t15 = qJD(4) * t17;
t20 = cos(qJ(2));
t24 = -t20 * pkin(1) - pkin(2);
t10 = -pkin(3) + t24;
t25 = t18 * t28;
t13 = t20 * t28;
t8 = t13 + qJD(3);
t26 = t10 * t16 + t17 * t25 + t19 * t8;
t23 = t19 * t25;
t11 = -0.2e1 * t25;
t4 = t17 * qJD(3) + (qJ(3) * t19 + t17 * t21) * qJD(4);
t3 = qJ(3) * t15 - t29;
t2 = -t23 + t17 * t8 + (t10 * t17 + t12 * t19) * qJD(4);
t1 = t12 * t15 - t26;
t5 = [0, 0, 0, 0, t11, -0.2e1 * t13, t11, 0.2e1 * t8, 0.2e1 * t12 * t8 + 0.2e1 * t24 * t25, 0, 0.2e1 * t2, -0.2e1 * t1; 0, 0, 0, 0, -t25, -t13, -t25, t22 + t13, -pkin(2) * t25 + t8 * qJ(3) + t12 * qJD(3), 0, -t23 + (qJD(3) + t8) * t17 + (t27 * t19 + (t10 + t21) * t17) * qJD(4), -t27 * t15 + t26 + t29; 0, 0, 0, 0, 0, 0, 0, t22, qJ(3) * t22, 0, 0.2e1 * t4, -0.2e1 * t3; 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t15, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
