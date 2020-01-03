% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:36
% EndTime: 2019-12-31 19:26:36
% DurationCPUTime: 0.17s
% Computational Cost: add. (102->40), mult. (278->64), div. (0->0), fcn. (177->6), ass. (0->32)
t25 = 2 * qJD(4);
t21 = sin(qJ(5));
t23 = cos(qJ(5));
t32 = qJD(5) * t23;
t19 = sin(pkin(8));
t20 = cos(pkin(8));
t24 = cos(qJ(2));
t34 = pkin(1) * qJD(2);
t30 = t24 * t34;
t22 = sin(qJ(2));
t31 = t22 * t34;
t8 = -t19 * t31 + t20 * t30;
t4 = qJD(4) + t8;
t16 = t24 * pkin(1) + pkin(2);
t36 = t20 * t22;
t26 = pkin(1) * t36 + t19 * t16;
t6 = qJ(4) + t26;
t37 = t4 * t21 + t6 * t32;
t15 = t19 * pkin(2) + qJ(4);
t35 = qJD(4) * t21 + t15 * t32;
t33 = qJD(5) * t21;
t29 = -t20 * pkin(2) - pkin(3);
t28 = -t19 * t22 * pkin(1) + t20 * t16;
t27 = -pkin(3) - t28;
t18 = qJD(4) * t23;
t14 = -pkin(7) + t29;
t12 = -0.2e1 * t21 * t32;
t9 = 0.2e1 * (t21 ^ 2 - t23 ^ 2) * qJD(5);
t7 = (t19 * t24 + t36) * t34;
t5 = -pkin(7) + t27;
t3 = t4 * t23;
t1 = [0, 0, 0, 0, -0.2e1 * t31, -0.2e1 * t30, 0.2e1 * t26 * t8 - 0.2e1 * t28 * t7, 0.2e1 * t7, 0.2e1 * t4, 0.2e1 * t27 * t7 + 0.2e1 * t6 * t4, t12, t9, 0, 0, 0, 0.2e1 * t37, -0.2e1 * t6 * t33 + 0.2e1 * t3; 0, 0, 0, 0, -t31, -t30, (t19 * t8 - t20 * t7) * pkin(2), t7, t25 + t8, t6 * qJD(4) + t4 * t15 + t7 * t29, t12, t9, 0, 0, 0, t35 + t37, t18 + t3 + (-t15 - t6) * t33; 0, 0, 0, 0, 0, 0, 0, 0, t25, t15 * t25, t12, t9, 0, 0, 0, 0.2e1 * t35, -0.2e1 * t15 * t33 + 0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t32, 0, t23 * t7 - t5 * t33, -t21 * t7 - t5 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t32, 0, -t14 * t33, -t14 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
