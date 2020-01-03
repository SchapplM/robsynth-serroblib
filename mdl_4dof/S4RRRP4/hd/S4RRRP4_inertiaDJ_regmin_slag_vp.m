% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S4RRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_inertiaDJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_inertiaDJ_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_inertiaDJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:41
% EndTime: 2019-12-31 17:15:42
% DurationCPUTime: 0.20s
% Computational Cost: add. (237->47), mult. (592->97), div. (0->0), fcn. (478->4), ass. (0->36)
t40 = qJD(2) + qJD(3);
t23 = sin(qJ(2));
t39 = pkin(5) + pkin(6);
t29 = qJD(2) * t39;
t13 = t23 * t29;
t25 = cos(qJ(2));
t14 = t25 * t29;
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t15 = t39 * t23;
t16 = t39 * t25;
t27 = -t24 * t15 - t22 * t16;
t3 = -t27 * qJD(3) + t24 * t13 + t22 * t14;
t38 = t24 * pkin(2);
t37 = t22 * t23;
t36 = qJD(2) * t23;
t35 = qJD(2) * t25;
t34 = qJD(3) * t24;
t33 = -0.2e1 * pkin(1) * qJD(2);
t32 = pkin(2) * t36;
t31 = qJD(3) * t22 * pkin(2);
t30 = pkin(2) * t34;
t21 = -t25 * pkin(2) - pkin(1);
t26 = t22 * t15 - t24 * t16;
t12 = t22 * t25 + t24 * t23;
t4 = t26 * qJD(3) + t22 * t13 - t24 * t14;
t20 = pkin(3) + t38;
t11 = -t24 * t25 + t37;
t9 = t40 * t12;
t8 = -t24 * t35 - t25 * t34 + t40 * t37;
t7 = t9 * pkin(3) + t32;
t6 = -t11 * qJ(4) - t26;
t5 = -t12 * qJ(4) + t27;
t2 = t8 * qJ(4) - t12 * qJD(4) + t4;
t1 = -t9 * qJ(4) - t11 * qJD(4) - t3;
t10 = [0, 0, 0, 0.2e1 * t23 * t35, 0.2e1 * (-t23 ^ 2 + t25 ^ 2) * qJD(2), 0, 0, 0, t23 * t33, t25 * t33, -0.2e1 * t12 * t8, 0.2e1 * t8 * t11 - 0.2e1 * t12 * t9, 0, 0, 0, 0.2e1 * t11 * t32 + 0.2e1 * t21 * t9, 0.2e1 * t12 * t32 - 0.2e1 * t21 * t8, -0.2e1 * t1 * t11 - 0.2e1 * t2 * t12 + 0.2e1 * t5 * t8 - 0.2e1 * t6 * t9, 0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * (t11 * pkin(3) + t21) * t7; 0, 0, 0, 0, 0, t35, -t36, 0, -pkin(5) * t35, pkin(5) * t36, 0, 0, -t8, -t9, 0, t4, t3, t20 * t8 + (-t22 * t9 + (-t11 * t24 + t12 * t22) * qJD(3)) * pkin(2), t2 * t20 + (t1 * t22 + (-t22 * t5 + t24 * t6) * qJD(3)) * pkin(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t31, -0.2e1 * t30, 0, 0.2e1 * (-t20 + t38) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, 0, t4, t3, pkin(3) * t8, t2 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t30, 0, -pkin(3) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
