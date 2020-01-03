% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:43
% EndTime: 2019-12-31 19:49:45
% DurationCPUTime: 0.54s
% Computational Cost: add. (257->61), mult. (752->100), div. (0->0), fcn. (488->6), ass. (0->52)
t40 = sin(pkin(8));
t41 = cos(pkin(8));
t45 = cos(qJ(2));
t57 = pkin(1) * qJD(2);
t55 = t45 * t57;
t43 = sin(qJ(2));
t56 = t43 * t57;
t19 = -t40 * t56 + t41 * t55;
t42 = sin(qJ(4));
t38 = t42 ^ 2;
t44 = cos(qJ(4));
t39 = t44 ^ 2;
t70 = (t38 + t39) * t19;
t68 = 0.2e1 * (-t38 + t39) * qJD(4);
t67 = 2 * qJD(5);
t34 = t45 * pkin(1) + pkin(2);
t60 = t41 * t43;
t58 = pkin(1) * t60 + t40 * t34;
t17 = pkin(7) + t58;
t66 = t70 * t17;
t65 = t41 * pkin(2);
t32 = t40 * pkin(2) + pkin(7);
t64 = t70 * t32;
t36 = t42 * qJD(4);
t37 = t44 * qJD(4);
t21 = -pkin(4) * t36 + qJ(5) * t37 + t42 * qJD(5);
t18 = (t40 * t45 + t60) * t57;
t7 = t18 - t21;
t63 = t21 - t7;
t48 = -t40 * t43 * pkin(1) + t41 * t34;
t16 = -pkin(3) - t48;
t59 = t16 * t37 + t18 * t42;
t33 = -pkin(3) - t65;
t54 = t33 * t36;
t53 = t33 * t37;
t52 = t42 * t37;
t51 = t32 * t36;
t50 = t32 * t37;
t49 = t16 * t36 - t18 * t44;
t47 = -t44 * pkin(4) - t42 * qJ(5);
t46 = -pkin(3) + t47;
t20 = t47 * qJD(4) + t44 * qJD(5);
t29 = -0.2e1 * t52;
t28 = 0.2e1 * t52;
t22 = t46 - t65;
t15 = t22 * t36;
t11 = t46 - t48;
t8 = t11 * t36;
t3 = t17 * t37 + t42 * t19;
t2 = t17 * t36 - t44 * t19;
t1 = 0.2e1 * t70;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t56, -0.2e1 * t55, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t18, -0.2e1 * t19, 0, -0.2e1 * t48 * t18 + 0.2e1 * t58 * t19, t28, t68, 0, t29, 0, 0, 0.2e1 * t49, 0.2e1 * t59, t1, 0.2e1 * t16 * t18 + 0.2e1 * t66, t28, 0, -t68, 0, 0, t29, -0.2e1 * t7 * t44 + 0.2e1 * t8, t1, -0.2e1 * t11 * t37 - 0.2e1 * t7 * t42, 0.2e1 * t11 * t7 + 0.2e1 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t55, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t19, 0, (-t18 * t41 + t19 * t40) * pkin(2), t28, t68, 0, t29, 0, 0, t49 + t54, t53 + t59, t70, t18 * t33 + t64, t28, 0, -t68, 0, 0, t29, t63 * t44 + t15 + t8, t70, t63 * t42 + (-t11 - t22) * t37, -t11 * t21 + t7 * t22 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t68, 0, t29, 0, 0, 0.2e1 * t54, 0.2e1 * t53, 0, 0, t28, 0, -t68, 0, 0, t29, 0.2e1 * t21 * t44 + 0.2e1 * t15, 0, 0.2e1 * t21 * t42 - 0.2e1 * t22 * t37, -0.2e1 * t22 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t36, 0, -t3, t2, 0, 0, 0, t37, 0, 0, t36, 0, -t3, t20, -t2, (-pkin(4) * t42 + qJ(5) * t44) * t19 + t20 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t36, 0, -t50, t51, 0, 0, 0, t37, 0, 0, t36, 0, -t50, t20, -t51, t20 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, t37, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, qJ(5) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t4;
