% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:33
% EndTime: 2019-12-31 17:03:33
% DurationCPUTime: 0.18s
% Computational Cost: add. (67->34), mult. (209->55), div. (0->0), fcn. (100->4), ass. (0->31)
t23 = 2 * qJD(3);
t19 = sin(qJ(2));
t10 = t19 * pkin(1) + qJ(3);
t18 = sin(qJ(4));
t20 = cos(qJ(4));
t29 = t20 * qJD(4);
t21 = cos(qJ(2));
t31 = pkin(1) * qJD(2);
t27 = t21 * t31;
t8 = qJD(3) + t27;
t34 = t10 * t29 + t8 * t18;
t33 = t10 * t8;
t28 = qJ(3) * qJD(4);
t32 = qJD(3) * t18 + t20 * t28;
t30 = t18 * qJD(4);
t12 = t19 * t31;
t26 = t18 * t29;
t25 = -t21 * pkin(1) - pkin(2);
t24 = t8 * qJ(3) + t10 * qJD(3);
t16 = t18 ^ 2;
t17 = t20 ^ 2;
t1 = (t16 + t17) * t12;
t22 = -pkin(2) - pkin(6);
t15 = qJ(3) * t23;
t14 = qJD(3) * t20;
t9 = -pkin(6) + t25;
t7 = -0.2e1 * t26;
t6 = 0.2e1 * t26;
t4 = t8 * t20;
t2 = 0.2e1 * (t16 - t17) * qJD(4);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t12, -0.2e1 * t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t12, 0.2e1 * t8, 0.2e1 * t25 * t12 + 0.2e1 * t33, t7, t2, 0, t6, 0, 0, 0.2e1 * t34, -0.2e1 * t10 * t30 + 0.2e1 * t4, -0.2e1 * t1, 0.2e1 * t9 * t1 + 0.2e1 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t23 + t27, -pkin(2) * t12 + t24, t7, t2, 0, t6, 0, 0, t32 + t34, t14 + t4 + (-qJ(3) - t10) * t30, -t1, t22 * t1 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t15, t7, t2, 0, t6, 0, 0, 0.2e1 * t32, -0.2e1 * t18 * t28 + 0.2e1 * t14, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, -t29, 0, t20 * t12 - t9 * t30, -t18 * t12 - t9 * t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, -t29, 0, -t22 * t30, -t22 * t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
