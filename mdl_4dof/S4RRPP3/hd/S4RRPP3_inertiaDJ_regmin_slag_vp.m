% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x16]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:51
% EndTime: 2019-12-31 16:57:51
% DurationCPUTime: 0.13s
% Computational Cost: add. (188->37), mult. (465->82), div. (0->0), fcn. (355->4), ass. (0->31)
t37 = 2 * qJD(4);
t36 = -qJ(3) - pkin(5);
t25 = sin(qJ(2));
t35 = t25 * qJD(2);
t26 = cos(qJ(2));
t34 = t26 * qJD(2);
t33 = -0.2e1 * pkin(1) * qJD(2);
t22 = pkin(2) * t35;
t29 = qJD(2) * t36;
t11 = t26 * qJD(3) + t25 * t29;
t23 = sin(pkin(6));
t24 = cos(pkin(6));
t28 = -t25 * qJD(3) + t26 * t29;
t5 = t23 * t11 - t24 * t28;
t6 = t24 * t11 + t23 * t28;
t17 = t36 * t26;
t30 = t36 * t25;
t8 = -t23 * t17 - t24 * t30;
t9 = -t24 * t17 + t23 * t30;
t32 = t8 * t5 + t9 * t6;
t31 = -t26 * pkin(2) - pkin(1);
t15 = t23 * t26 + t24 * t25;
t12 = t15 * qJD(2);
t13 = -t23 * t35 + t24 * t34;
t14 = t23 * t25 - t24 * t26;
t27 = -0.2e1 * t9 * t12 + 0.2e1 * t8 * t13 - 0.2e1 * t6 * t14 + 0.2e1 * t5 * t15;
t21 = -t24 * pkin(2) - pkin(3);
t19 = t23 * pkin(2) + qJ(4);
t7 = t14 * pkin(3) - t15 * qJ(4) + t31;
t2 = t12 * pkin(3) - t13 * qJ(4) - t15 * qJD(4) + t22;
t1 = [0, 0, 0, 0.2e1 * t25 * t34, 0.2e1 * (-t25 ^ 2 + t26 ^ 2) * qJD(2), 0, 0, 0, t25 * t33, t26 * t33, t27, 0.2e1 * t31 * t22 + 0.2e1 * t32, 0.2e1 * t7 * t12 + 0.2e1 * t2 * t14, t27, -0.2e1 * t7 * t13 - 0.2e1 * t2 * t15, 0.2e1 * t7 * t2 + 0.2e1 * t32; 0, 0, 0, 0, 0, t34, -t35, 0, -pkin(5) * t34, pkin(5) * t35, (-t12 * t23 - t13 * t24) * pkin(2), (t23 * t6 - t24 * t5) * pkin(2), -t5, -qJD(4) * t14 - t19 * t12 + t21 * t13, t6, t9 * qJD(4) + t6 * t19 + t5 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t19 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t12, 0, -t13, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
