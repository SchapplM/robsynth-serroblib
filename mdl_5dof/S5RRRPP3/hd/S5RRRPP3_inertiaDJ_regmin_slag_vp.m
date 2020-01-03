% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:41
% EndTime: 2019-12-31 20:53:43
% DurationCPUTime: 0.42s
% Computational Cost: add. (356->97), mult. (823->154), div. (0->0), fcn. (492->4), ass. (0->66)
t53 = cos(qJ(3));
t44 = t53 * qJD(4);
t51 = sin(qJ(3));
t72 = t51 * qJ(4);
t59 = -t53 * pkin(3) - t72;
t80 = t59 * qJD(3) + t44;
t50 = -pkin(3) - qJ(5);
t79 = t50 * t53 - t72;
t71 = t51 * qJD(3);
t61 = pkin(3) * t71 - t51 * qJD(4);
t69 = qJD(3) * qJ(4);
t4 = qJ(5) * t71 + (-qJD(5) - t69) * t53 + t61;
t52 = sin(qJ(2));
t73 = pkin(1) * qJD(2);
t43 = t52 * t73;
t1 = t43 + t4;
t78 = -t1 - t4;
t54 = cos(qJ(2));
t77 = t54 * pkin(1);
t39 = t52 * pkin(1) + pkin(7);
t76 = pkin(4) + t39;
t40 = -pkin(2) - t77;
t45 = t53 * qJD(3);
t74 = t40 * t45 + t51 * t43;
t70 = qJ(4) * qJD(4);
t46 = t51 * pkin(4);
t24 = t51 * t39 + t46;
t64 = t54 * t73;
t35 = t53 * t64;
t6 = -t76 * t71 + t35;
t60 = t51 * t64;
t7 = t76 * t45 + t60;
t68 = t24 * t45 + t7 * t51 + t6 * t53;
t28 = (-pkin(4) - pkin(7)) * t71;
t41 = pkin(7) * t45;
t29 = pkin(4) * t45 + t41;
t33 = t51 * pkin(7) + t46;
t67 = t28 * t53 + t29 * t51 + t33 * t45;
t66 = pkin(2) * t71;
t65 = pkin(2) * t45;
t63 = pkin(7) * t71;
t32 = -pkin(2) + t59;
t23 = t32 - t77;
t62 = qJD(3) * (-t23 - t32);
t48 = t51 ^ 2;
t49 = t53 ^ 2;
t17 = (t48 + t49) * t64;
t58 = t40 * t71 - t53 * t43;
t22 = -pkin(2) + t79;
t19 = -t53 * t69 + t61;
t55 = 0.2e1 * qJD(4);
t47 = t53 * pkin(4);
t37 = 0.2e1 * t51 * t45;
t34 = t53 * pkin(7) + t47;
t27 = 0.2e1 * (-t48 + t49) * qJD(3);
t25 = t53 * t39 + t47;
t15 = t22 - t77;
t14 = t39 * t45 + t60;
t13 = t39 * t71 - t35;
t12 = t22 * t71;
t11 = t19 * t53;
t10 = t19 + t43;
t9 = t79 * qJD(3) - qJD(5) * t51 + t44;
t8 = t15 * t71;
t5 = t10 * t53;
t2 = [0, 0, 0, 0, -0.2e1 * t43, -0.2e1 * t64, t37, t27, 0, 0, 0, 0.2e1 * t58, 0.2e1 * t74, 0.2e1 * t17, -0.2e1 * t23 * t71 + 0.2e1 * t5, -0.2e1 * t10 * t51 - 0.2e1 * t23 * t45, 0.2e1 * t23 * t10 + 0.2e1 * t39 * t17, -0.2e1 * t25 * t71 + 0.2e1 * t68, -0.2e1 * t1 * t51 - 0.2e1 * t15 * t45, -0.2e1 * t1 * t53 + 0.2e1 * t8, 0.2e1 * t15 * t1 + 0.2e1 * t24 * t7 + 0.2e1 * t25 * t6; 0, 0, 0, 0, -t43, -t64, t37, t27, 0, 0, 0, t58 - t66, -t65 + t74, t17, t51 * t62 + t11 + t5, (-t10 - t19) * t51 + t53 * t62, pkin(7) * t17 + t10 * t32 + t23 * t19, (-t25 - t34) * t71 + t67 + t68, t78 * t51 + (-t15 - t22) * t45, t78 * t53 + t12 + t8, t1 * t22 + t15 * t4 + t24 * t29 + t25 * t28 + t7 * t33 + t6 * t34; 0, 0, 0, 0, 0, 0, t37, t27, 0, 0, 0, -0.2e1 * t66, -0.2e1 * t65, 0, -0.2e1 * t32 * t71 + 0.2e1 * t11, -0.2e1 * t19 * t51 - 0.2e1 * t32 * t45, 0.2e1 * t32 * t19, -0.2e1 * t34 * t71 + 0.2e1 * t67, -0.2e1 * t22 * t45 - 0.2e1 * t4 * t51, -0.2e1 * t4 * t53 + 0.2e1 * t12, 0.2e1 * t22 * t4 + 0.2e1 * t34 * t28 + 0.2e1 * t33 * t29; 0, 0, 0, 0, 0, 0, 0, 0, t45, -t71, 0, -t14, t13, t80, t14, -t13, (-pkin(3) * t51 + qJ(4) * t53) * t64 + t80 * t39, t9, t6, -t7, t6 * qJ(4) + t25 * qJD(4) - t24 * qJD(5) + t7 * t50; 0, 0, 0, 0, 0, 0, 0, 0, t45, -t71, 0, -t41, t63, t80, t41, -t63, t80 * pkin(7), t9, t28, -t29, t28 * qJ(4) + t34 * qJD(4) - t33 * qJD(5) + t29 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0.2e1 * t70, 0, t55, 0.2e1 * qJD(5), -0.2e1 * t50 * qJD(5) + 0.2e1 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, t14, t45, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, t41, t45, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t2;
