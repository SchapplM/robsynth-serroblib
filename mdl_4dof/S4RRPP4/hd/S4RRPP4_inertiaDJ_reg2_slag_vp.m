% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPP4
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
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPP4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_inertiaDJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:18
% EndTime: 2019-12-31 16:59:19
% DurationCPUTime: 0.23s
% Computational Cost: add. (86->37), mult. (240->70), div. (0->0), fcn. (133->2), ass. (0->31)
t16 = cos(qJ(2));
t15 = sin(qJ(2));
t28 = t15 * qJ(3);
t31 = pkin(2) + pkin(3);
t32 = t31 * t16 + t28;
t6 = 0.2e1 * (-t15 ^ 2 + t16 ^ 2) * qJD(2);
t18 = 2 * qJD(3);
t30 = pkin(5) - qJ(4);
t29 = qJ(3) * t16;
t27 = t15 * qJD(2);
t26 = t15 * qJD(3);
t13 = t16 * qJD(2);
t25 = t16 * qJD(3);
t24 = -0.2e1 * pkin(1) * qJD(2);
t23 = pkin(5) * t27;
t22 = pkin(5) * t13;
t21 = t15 * t13;
t10 = t30 * t16;
t20 = -t16 * pkin(2) - t28;
t19 = t20 * qJD(2) + t25;
t14 = qJ(3) * t18;
t12 = -0.2e1 * t21;
t11 = 0.2e1 * t21;
t9 = t30 * t15;
t7 = -pkin(1) + t20;
t5 = pkin(1) + t32;
t4 = -t26 + (pkin(2) * t15 - t29) * qJD(2);
t3 = qJD(2) * t10 - t15 * qJD(4);
t2 = -t16 * qJD(4) - t30 * t27;
t1 = t26 + (-t31 * t15 + t29) * qJD(2);
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t6, 0, t12, 0, 0, t15 * t24, t16 * t24, 0, 0, t11, 0, -t6, 0, 0, t12, -0.2e1 * t4 * t16 + 0.2e1 * t7 * t27, 0, -0.2e1 * t7 * t13 - 0.2e1 * t4 * t15, 0.2e1 * t7 * t4, t11, -t6, 0, t12, 0, 0, 0.2e1 * t1 * t16 - 0.2e1 * t5 * t27, 0.2e1 * t1 * t15 + 0.2e1 * t5 * t13, -0.2e1 * t3 * t15 - 0.2e1 * t2 * t16 + 0.2e1 * (t10 * t15 - t16 * t9) * qJD(2), 0.2e1 * t5 * t1 + 0.2e1 * t10 * t2 + 0.2e1 * t9 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, -t27, 0, -t22, t23, 0, 0, 0, t13, 0, 0, t27, 0, -t22, t19, -t23, t19 * pkin(5), 0, 0, -t13, 0, -t27, 0, -t3, t2, t32 * qJD(2) - t25, t2 * qJ(3) + t10 * qJD(3) - t3 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t14, 0, 0, 0, 0, 0, 0, 0, t18, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t13, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
