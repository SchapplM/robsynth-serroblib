% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S4RRPR4
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S4RRPR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_inertiaDJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_inertiaDJ_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:34
% DurationCPUTime: 0.34s
% Computational Cost: add. (336->57), mult. (816->102), div. (0->0), fcn. (658->6), ass. (0->50)
t50 = sin(pkin(7));
t51 = cos(pkin(7));
t65 = t50 ^ 2 + t51 ^ 2;
t76 = t65 * qJD(3);
t77 = 0.2e1 * t76;
t54 = cos(qJ(2));
t64 = pkin(1) * qJD(2);
t63 = t54 * t64;
t38 = qJD(3) + t63;
t75 = t65 * t38;
t74 = t54 * pkin(1);
t73 = cos(qJ(4));
t58 = qJD(4) * t73;
t52 = sin(qJ(4));
t70 = t52 * t50;
t25 = qJD(4) * t70 - t51 * t58;
t43 = -t51 * pkin(3) - pkin(2);
t72 = t43 * t25;
t32 = t73 * t50 + t52 * t51;
t26 = t32 * qJD(4);
t71 = t43 * t26;
t61 = t73 * t51;
t31 = -t61 + t70;
t35 = t43 - t74;
t53 = sin(qJ(2));
t44 = t53 * t64;
t69 = t35 * t26 + t31 * t44;
t68 = -t35 * t25 + t32 * t44;
t42 = t53 * pkin(1) + qJ(3);
t27 = (-pkin(6) - t42) * t50;
t47 = t51 * pkin(6);
t28 = t51 * t42 + t47;
t14 = t73 * t27 - t52 * t28;
t15 = t52 * t27 + t73 * t28;
t3 = -t27 * t58 - t38 * t61 + (qJD(4) * t28 + t38 * t50) * t52;
t4 = -t15 * qJD(4) - t32 * t38;
t62 = t14 * t25 - t15 * t26 + t3 * t31 - t4 * t32;
t36 = (-pkin(6) - qJ(3)) * t50;
t37 = t51 * qJ(3) + t47;
t10 = -t36 * t58 - qJD(3) * t61 + (qJD(3) * t50 + qJD(4) * t37) * t52;
t19 = t52 * t36 + t73 * t37;
t11 = -t32 * qJD(3) - t19 * qJD(4);
t18 = t73 * t36 - t52 * t37;
t59 = t10 * t31 - t11 * t32 + t18 * t25 - t19 * t26;
t56 = t50 * t44;
t55 = t51 * t44;
t17 = -0.2e1 * t32 * t25;
t16 = 0.2e1 * t31 * t26;
t5 = 0.2e1 * t25 * t31 - 0.2e1 * t32 * t26;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t44, -0.2e1 * t63, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t55, 0.2e1 * t56, 0.2e1 * t75, 0.2e1 * (-pkin(2) - t74) * t44 + 0.2e1 * t42 * t75, t17, t5, 0, t16, 0, 0, 0.2e1 * t69, 0.2e1 * t68, 0.2e1 * t62, 0.2e1 * t14 * t4 - 0.2e1 * t15 * t3 + 0.2e1 * t35 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t63, 0, 0, 0, 0, 0, 0, 0, 0, -t55, t56, t76 + t75, -pkin(2) * t44 + qJ(3) * t75 + t42 * t76, t17, t5, 0, t16, 0, 0, t69 + t71, t68 - t72, t59 + t62, -t15 * t10 + t14 * t11 + t4 * t18 - t3 * t19 + t43 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, qJ(3) * t77, t17, t5, 0, t16, 0, 0, 0.2e1 * t71, -0.2e1 * t72, 0.2e1 * t59, -0.2e1 * t19 * t10 + 0.2e1 * t18 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, 0, t26, -t25, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, -t26, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, -t26, 0, t11, t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
