% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:34
% EndTime: 2019-12-31 17:01:34
% DurationCPUTime: 0.19s
% Computational Cost: add. (90->30), mult. (324->60), div. (0->0), fcn. (205->6), ass. (0->33)
t25 = cos(qJ(4));
t18 = t25 * qJD(4);
t23 = sin(qJ(4));
t26 = cos(qJ(2));
t17 = t26 * pkin(1) + pkin(2);
t21 = sin(pkin(7));
t22 = cos(pkin(7));
t24 = sin(qJ(2));
t27 = -t21 * t24 * pkin(1) + t22 * t17;
t5 = -pkin(3) - t27;
t35 = pkin(1) * qJD(2);
t37 = t22 * t24;
t7 = (t21 * t26 + t37) * t35;
t38 = t5 * t18 + t7 * t23;
t36 = pkin(1) * t37 + t21 * t17;
t34 = t23 * qJD(4);
t33 = t24 * t35;
t32 = t26 * t35;
t16 = -t22 * pkin(2) - pkin(3);
t31 = t16 * t34;
t30 = t16 * t18;
t29 = t23 * t18;
t28 = -t7 * t25 + t5 * t34;
t19 = t23 ^ 2;
t20 = t25 ^ 2;
t8 = -t21 * t33 + t22 * t32;
t1 = (t19 + t20) * t8;
t15 = t21 * pkin(2) + pkin(6);
t13 = -0.2e1 * t29;
t12 = 0.2e1 * t29;
t9 = 0.2e1 * (-t19 + t20) * qJD(4);
t6 = pkin(6) + t36;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t33, -0.2e1 * t32, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t7, -0.2e1 * t8, 0, -0.2e1 * t27 * t7 + 0.2e1 * t36 * t8, t12, t9, 0, t13, 0, 0, 0.2e1 * t28, 0.2e1 * t38, 0.2e1 * t1, 0.2e1 * t6 * t1 + 0.2e1 * t5 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t32, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, (t21 * t8 - t22 * t7) * pkin(2), t12, t9, 0, t13, 0, 0, t28 + t31, t30 + t38, t1, t15 * t1 + t7 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t9, 0, t13, 0, 0, 0.2e1 * t31, 0.2e1 * t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t34, 0, -t6 * t18 - t23 * t8, -t25 * t8 + t6 * t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t34, 0, -t15 * t18, t15 * t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
