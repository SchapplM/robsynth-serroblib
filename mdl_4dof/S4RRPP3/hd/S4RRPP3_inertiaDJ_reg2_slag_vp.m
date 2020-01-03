% Calculate inertial parameters regressor of joint inertia matrix time derivative for
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
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:51
% EndTime: 2019-12-31 16:57:52
% DurationCPUTime: 0.24s
% Computational Cost: add. (230->43), mult. (606->97), div. (0->0), fcn. (481->4), ass. (0->35)
t41 = 2 * qJD(4);
t40 = -qJ(3) - pkin(5);
t27 = sin(qJ(2));
t39 = t27 * qJD(2);
t28 = cos(qJ(2));
t38 = t28 * qJD(2);
t25 = sin(pkin(6));
t26 = cos(pkin(6));
t16 = t25 * t28 + t26 * t27;
t13 = t16 * qJD(2);
t15 = t25 * t27 - t26 * t28;
t37 = 0.2e1 * t15 * t13;
t36 = -0.2e1 * pkin(1) * qJD(2);
t24 = pkin(2) * t39;
t18 = t40 * t28;
t33 = t40 * t27;
t10 = -t26 * t18 + t25 * t33;
t32 = qJD(2) * t40;
t12 = t28 * qJD(3) + t27 * t32;
t30 = -t27 * qJD(3) + t28 * t32;
t5 = t25 * t12 - t26 * t30;
t6 = t26 * t12 + t25 * t30;
t9 = -t25 * t18 - t26 * t33;
t35 = t10 * t6 + t9 * t5;
t34 = t27 * t38;
t23 = -t28 * pkin(2) - pkin(1);
t14 = -t25 * t39 + t26 * t38;
t31 = t16 * t13 + t14 * t15;
t29 = -0.2e1 * t10 * t13 + 0.2e1 * t9 * t14 - 0.2e1 * t6 * t15 + 0.2e1 * t5 * t16;
t22 = -t26 * pkin(2) - pkin(3);
t20 = t25 * pkin(2) + qJ(4);
t8 = 0.2e1 * t16 * t14;
t7 = t15 * pkin(3) - t16 * qJ(4) + t23;
t2 = t13 * pkin(3) - t14 * qJ(4) - t16 * qJD(4) + t24;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t34, 0.2e1 * (-t27 ^ 2 + t28 ^ 2) * qJD(2), 0, -0.2e1 * t34, 0, 0, t27 * t36, t28 * t36, 0, 0, t8, -0.2e1 * t31, 0, t37, 0, 0, 0.2e1 * t23 * t13 + 0.2e1 * t15 * t24, 0.2e1 * t23 * t14 + 0.2e1 * t16 * t24, t29, 0.2e1 * t23 * t24 + 0.2e1 * t35, t8, 0, 0.2e1 * t31, 0, 0, t37, 0.2e1 * t7 * t13 + 0.2e1 * t2 * t15, t29, -0.2e1 * t7 * t14 - 0.2e1 * t2 * t16, 0.2e1 * t7 * t2 + 0.2e1 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, -t39, 0, -pkin(5) * t38, pkin(5) * t39, 0, 0, 0, 0, t14, 0, -t13, 0, -t5, -t6, (-t13 * t25 - t14 * t26) * pkin(2), (t25 * t6 - t26 * t5) * pkin(2), 0, t14, 0, 0, t13, 0, -t5, -qJD(4) * t15 - t20 * t13 + t22 * t14, t6, t10 * qJD(4) + t6 * t20 + t5 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t20 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, t24, 0, 0, 0, 0, 0, 0, t13, 0, -t14, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
