% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RPRR2
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
% MMD_reg [((4+1)*4/2)x14]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RPRR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:13
% EndTime: 2019-12-31 16:48:13
% DurationCPUTime: 0.13s
% Computational Cost: add. (54->18), mult. (157->33), div. (0->0), fcn. (100->6), ass. (0->21)
t11 = cos(pkin(7)) * pkin(1) + pkin(2);
t15 = sin(qJ(3));
t17 = cos(qJ(3));
t25 = pkin(1) * sin(pkin(7));
t27 = -t17 * t11 + t15 * t25;
t16 = cos(qJ(4));
t12 = t16 * qJD(4);
t14 = sin(qJ(4));
t18 = t15 * t11 + t17 * t25;
t5 = t18 * qJD(3);
t6 = -pkin(3) + t27;
t26 = t6 * t12 + t5 * t14;
t23 = t14 * qJD(4);
t21 = pkin(3) * t23;
t20 = pkin(3) * t12;
t19 = -t5 * t16 + t6 * t23;
t10 = 0.2e1 * t14 * t12;
t8 = 0.2e1 * (-t14 ^ 2 + t16 ^ 2) * qJD(4);
t7 = pkin(6) + t18;
t4 = t27 * qJD(3);
t1 = [0, 0, 0, 0, 0, -0.2e1 * t5, 0.2e1 * t4, t10, t8, 0, 0, 0, 0.2e1 * t19, 0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t5, t4, t10, t8, 0, 0, 0, t19 - t21, -t20 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t10, t8, 0, 0, 0, -0.2e1 * t21, -0.2e1 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t23, 0, -t7 * t12 + t14 * t4, t16 * t4 + t7 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t23, 0, -pkin(6) * t12, pkin(6) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
