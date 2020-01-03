% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4PRRP5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaDJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:22
% EndTime: 2019-12-31 16:29:23
% DurationCPUTime: 0.23s
% Computational Cost: add. (69->29), mult. (252->66), div. (0->0), fcn. (164->4), ass. (0->34)
t34 = 2 * qJD(3);
t21 = cos(qJ(3));
t33 = t21 * pkin(3);
t32 = -qJ(4) - pkin(5);
t19 = sin(qJ(3));
t17 = t19 ^ 2;
t18 = t21 ^ 2;
t31 = t17 + t18;
t30 = t19 * qJD(3);
t20 = sin(qJ(2));
t29 = t20 * qJD(2);
t16 = t21 * qJD(3);
t22 = cos(qJ(2));
t28 = t22 * qJD(2);
t27 = -2 * pkin(2) * qJD(3);
t26 = pkin(3) * t30;
t25 = t19 * t16;
t24 = t31 * t22;
t6 = -t20 * t16 - t19 * t28;
t10 = t32 * t19;
t11 = t32 * t21;
t2 = -t21 * qJD(4) - t32 * t30;
t3 = -t19 * qJD(4) + t32 * t16;
t23 = -t3 * t19 - t2 * t21 + (-t10 * t21 + t11 * t19) * qJD(3);
t14 = -pkin(2) - t33;
t13 = -0.2e1 * t25;
t12 = 0.2e1 * t25;
t9 = (-t17 + t18) * t34;
t8 = -t21 * t29 - t22 * t30;
t7 = -t22 * t16 + t19 * t29;
t5 = qJD(2) * t24;
t4 = t20 * t30 - t21 * t28;
t1 = 0.2e1 * (-0.1e1 + t31) * t20 * t28;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t28, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, t5, (-pkin(2) * t20 + pkin(5) * t24) * qJD(2), 0, 0, 0, 0, 0, 0, t8, t7, t5, (-t26 + (-t10 * t19 - t11 * t21) * qJD(2)) * t22 + (qJD(2) * t14 + t23) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t9, 0, t13, 0, 0, t19 * t27, t21 * t27, 0, 0, t12, t9, 0, t13, 0, 0, 0.2e1 * (t14 - t33) * t30, (pkin(3) * t17 + t14 * t21) * t34, 0.2e1 * t23, 0.2e1 * t10 * t3 + 0.2e1 * t11 * t2 + 0.2e1 * t14 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t4, 0, 0, 0, 0, 0, 0, 0, 0, t6, t4, 0, t6 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t30, 0, -pkin(5) * t16, pkin(5) * t30, 0, 0, 0, 0, t16, 0, -t30, 0, t3, t2, -pkin(3) * t16, t3 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t16, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
