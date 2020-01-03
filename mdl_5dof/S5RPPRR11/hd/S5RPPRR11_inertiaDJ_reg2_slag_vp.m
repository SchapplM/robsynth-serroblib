% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR11_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:54
% EndTime: 2019-12-31 18:05:56
% DurationCPUTime: 0.62s
% Computational Cost: add. (384->76), mult. (791->161), div. (0->0), fcn. (539->4), ass. (0->65)
t22 = -pkin(6) + qJ(2);
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t40 = pkin(4) * t27 + pkin(7) * t25;
t63 = qJD(5) * t25;
t73 = t40 * qJD(4) - t22 * t63 + qJD(3);
t23 = pkin(1) + qJ(3);
t70 = t25 * pkin(4);
t39 = -t27 * pkin(7) + t70;
t34 = t39 + t23;
t57 = t27 * qJD(4);
t46 = t22 * t57;
t36 = -t25 * qJD(2) - t46;
t72 = -qJD(5) * t34 + t36;
t24 = sin(qJ(5));
t18 = t24 ^ 2;
t26 = cos(qJ(5));
t20 = t26 ^ 2;
t68 = t18 - t20;
t44 = qJD(5) * t68;
t19 = t25 ^ 2;
t21 = t27 ^ 2;
t43 = (t19 - t21) * qJD(4);
t58 = t25 * qJD(4);
t35 = t27 * qJD(2) - t22 * t58;
t71 = t40 * qJD(5) - t35;
t69 = t22 * t25;
t67 = t18 + t20;
t65 = t19 + t21;
t64 = qJD(5) * t24;
t17 = qJD(5) * t26;
t62 = qJD(5) * t27;
t61 = t21 * qJD(2);
t60 = t23 * qJD(3);
t56 = qJ(2) * qJD(2);
t55 = -0.2e1 * qJD(5) * pkin(4);
t54 = t24 * t69;
t53 = t26 * t69;
t52 = t24 * t62;
t50 = t26 * t62;
t49 = t24 * t17;
t48 = t26 * t58;
t47 = t25 * t57;
t45 = t67 * t27;
t14 = 0.2e1 * t47;
t42 = t24 * t48;
t41 = t21 * t49;
t4 = t26 * t34 - t54;
t5 = t24 * t34 + t53;
t38 = t24 * t5 + t26 * t4;
t37 = t24 * t4 - t26 * t5;
t31 = t39 * qJD(4) - t22 * t62;
t1 = -t73 * t24 + t72 * t26;
t2 = t72 * t24 + t73 * t26;
t30 = t37 * qJD(5) + t1 * t24 - t2 * t26;
t29 = -t38 * qJD(5) - t1 * t26 - t2 * t24;
t28 = 0.2e1 * qJD(2);
t15 = t22 * t61;
t12 = t65 * qJD(2);
t11 = t24 * t58 - t50;
t10 = t25 * t17 + t24 * t57;
t9 = -t48 - t52;
t8 = t24 * t63 - t26 * t57;
t3 = t27 * t44 + t42;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0.2e1 * t56, 0, 0, 0, 0, 0, 0, 0, t28, 0.2e1 * qJD(3), 0.2e1 * t56 + 0.2e1 * t60, -0.2e1 * t47, 0.2e1 * t43, 0, t14, 0, 0, 0.2e1 * qJD(3) * t25 + 0.2e1 * t23 * t57, 0.2e1 * qJD(3) * t27 - 0.2e1 * t23 * t58, -0.2e1 * t12, 0.2e1 * t19 * t22 * qJD(2) + 0.2e1 * t15 + 0.2e1 * t60, -0.2e1 * t20 * t47 - 0.2e1 * t41, 0.2e1 * t21 * t44 + 0.4e1 * t27 * t42, -0.2e1 * t25 * t52 - 0.2e1 * t26 * t43, -0.2e1 * t18 * t47 + 0.2e1 * t41, 0.2e1 * t24 * t43 - 0.2e1 * t25 * t50, t14, 0.2e1 * t2 * t25 + 0.2e1 * (-t24 * qJD(2) - t22 * t17) * t21 + 0.2e1 * (t4 + 0.2e1 * t54) * t57, 0.2e1 * t1 * t25 + 0.2e1 * (-qJD(2) * t26 + t22 * t64) * t21 + 0.2e1 * (-t5 + 0.2e1 * t53) * t57, 0.2e1 * t30 * t27 + 0.2e1 * t38 * t58, -0.2e1 * t22 ^ 2 * t47 - 0.2e1 * t5 * t1 + 0.2e1 * t4 * t2 + 0.2e1 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), 0, 0, 0, 0, 0, 0, -t57, t58, 0, -qJD(3), 0, 0, 0, 0, 0, 0, t8, t10, -t67 * t58, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, -t65 * t17, t65 * t64, 0, t61 - t37 * t57 + (t29 - 0.2e1 * t46) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t67) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, 0, -t57, 0, t35, t36, 0, 0, -t3, -0.4e1 * t27 * t49 + t68 * t58, t10, t3, -t8, 0, t31 * t24 - t71 * t26, t71 * t24 + t31 * t26, t29, t35 * pkin(4) + t29 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t57, 0, 0, 0, 0, 0, 0, 0, 0, t9, t11, qJD(4) * t45, (pkin(7) * t45 - t70) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t49, -0.2e1 * t44, 0, -0.2e1 * t49, 0, 0, t24 * t55, t26 * t55, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, t11, t57, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t8, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t64, 0, -pkin(7) * t17, pkin(7) * t64, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
