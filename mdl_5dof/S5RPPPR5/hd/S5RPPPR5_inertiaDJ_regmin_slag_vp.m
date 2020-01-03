% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPPR5_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:30
% EndTime: 2019-12-31 17:46:31
% DurationCPUTime: 0.15s
% Computational Cost: add. (84->24), mult. (195->60), div. (0->0), fcn. (162->6), ass. (0->26)
t15 = sin(pkin(8));
t17 = cos(pkin(8));
t18 = cos(pkin(7));
t27 = t18 * qJD(2);
t9 = -qJD(4) + t27;
t25 = (t15 ^ 2 + t17 ^ 2) * t9;
t32 = 0.2e1 * qJD(2);
t16 = sin(pkin(7));
t21 = -pkin(1) - pkin(2);
t30 = t18 * qJ(2) + t16 * t21;
t8 = -qJ(4) + t30;
t31 = pkin(6) - t8;
t28 = t16 * qJD(2);
t26 = 0.2e1 * t28;
t24 = -t16 * qJ(2) + t18 * t21;
t23 = pkin(3) - t24;
t19 = sin(qJ(5));
t20 = cos(qJ(5));
t22 = t20 * t15 + t19 * t17;
t6 = t19 * t15 - t20 * t17;
t5 = t22 * qJD(5);
t4 = t6 * qJD(5);
t3 = t17 * pkin(4) + t23;
t2 = t31 * t17;
t1 = t31 * t15;
t7 = [0, 0, 0, 0, t32, qJ(2) * t32, t26, 0.2e1 * t27, (-t24 * t16 + t30 * t18) * t32, t17 * t26, -0.2e1 * t15 * t28, -0.2e1 * t25, 0.2e1 * t23 * t28 + 0.2e1 * t8 * t25, -0.2e1 * t22 * t4, -0.2e1 * t22 * t5 + 0.2e1 * t4 * t6, 0, 0, 0, -0.2e1 * t6 * t28 - 0.2e1 * t3 * t5, -0.2e1 * t22 * t28 + 0.2e1 * t3 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t25 - t27) * t16, 0, 0, 0, 0, 0, t18 * t5, -t18 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, -t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, -t22 * t9 + (-t1 * t19 + t2 * t20) * qJD(5), t6 * t9 + (-t1 * t20 - t19 * t2) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t4, t16 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t7;
