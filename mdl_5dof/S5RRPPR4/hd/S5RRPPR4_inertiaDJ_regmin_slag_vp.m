% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:48
% EndTime: 2019-12-31 19:27:49
% DurationCPUTime: 0.19s
% Computational Cost: add. (147->50), mult. (324->89), div. (0->0), fcn. (210->6), ass. (0->39)
t37 = 2 * qJD(3);
t35 = cos(qJ(2));
t44 = pkin(1) * qJD(2);
t27 = t35 * t44;
t19 = t27 + qJD(3);
t30 = sin(pkin(8));
t31 = cos(pkin(8));
t33 = sin(qJ(2));
t40 = t33 * t44;
t7 = t31 * t19 + t30 * t40;
t39 = -t35 * pkin(1) - pkin(2);
t23 = -pkin(3) + t39;
t25 = t33 * pkin(1) + qJ(3);
t5 = t30 * t23 + t31 * t25;
t36 = -pkin(2) - pkin(3);
t12 = t31 * qJ(3) + t30 * t36;
t32 = sin(qJ(5));
t29 = qJD(5) * t32;
t34 = cos(qJ(5));
t43 = qJD(5) * t34;
t42 = t30 * qJD(3);
t41 = t31 * qJD(3);
t4 = t31 * t23 - t30 * t25;
t2 = pkin(4) - t4;
t11 = -t30 * qJ(3) + t31 * t36;
t9 = pkin(4) - t11;
t38 = qJD(5) * (-t2 - t9);
t24 = -0.2e1 * t40;
t22 = t31 * t43;
t21 = t31 * t29;
t20 = t34 * t42;
t18 = 0.2e1 * t32 * t43;
t17 = t31 * t40;
t10 = -pkin(7) + t12;
t8 = 0.2e1 * (-t32 ^ 2 + t34 ^ 2) * qJD(5);
t6 = t30 * t19 - t17;
t3 = -pkin(7) + t5;
t1 = t6 * t34;
t13 = [0, 0, 0, 0, t24, -0.2e1 * t27, t24, 0.2e1 * t19, 0.2e1 * t25 * t19 + 0.2e1 * t39 * t40, 0.2e1 * t6, 0.2e1 * t7, -0.2e1 * t4 * t6 + 0.2e1 * t5 * t7, t18, t8, 0, 0, 0, -0.2e1 * t2 * t29 + 0.2e1 * t1, -0.2e1 * t2 * t43 - 0.2e1 * t6 * t32; 0, 0, 0, 0, -t40, -t27, -t40, t37 + t27, -pkin(2) * t40 + t19 * qJ(3) + t25 * qJD(3), -t17 + (qJD(3) + t19) * t30, t41 + t7, -t6 * t11 + t7 * t12 + (-t30 * t4 + t31 * t5) * qJD(3), t18, t8, 0, 0, 0, t32 * t38 + t1 + t20, (-t6 - t42) * t32 + t34 * t38; 0, 0, 0, 0, 0, 0, 0, t37, qJ(3) * t37, 0.2e1 * t42, 0.2e1 * t41, (-t11 * t30 + t12 * t31) * t37, t18, t8, 0, 0, 0, -0.2e1 * t9 * t29 + 0.2e1 * t20, -0.2e1 * t32 * t42 - 0.2e1 * t9 * t43; 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, t7 * t30 - t6 * t31, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t29, 0, -t3 * t43 - t32 * t7, t3 * t29 - t34 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t29, 0, -t10 * t43 - t32 * t41, t10 * t29 - t34 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t43, t30 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
