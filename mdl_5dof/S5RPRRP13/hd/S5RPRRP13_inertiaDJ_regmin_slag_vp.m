% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP13_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_inertiaDJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:37
% EndTime: 2019-12-31 18:59:40
% DurationCPUTime: 0.77s
% Computational Cost: add. (542->122), mult. (1179->215), div. (0->0), fcn. (814->4), ass. (0->84)
t33 = sin(qJ(3));
t35 = cos(qJ(3));
t88 = pkin(7) * t35;
t49 = pkin(3) * t33 - t88;
t21 = qJ(2) + t49;
t32 = sin(qJ(4));
t19 = t32 * t21;
t34 = cos(qJ(4));
t36 = -pkin(1) - pkin(6);
t24 = t34 * t33 * t36;
t42 = t24 + t19;
t28 = t32 ^ 2;
t30 = t34 ^ 2;
t83 = t28 - t30;
t55 = t83 * qJD(4);
t44 = pkin(4) * t32 - qJ(5) * t34;
t41 = -t36 + t44;
t8 = t41 * t35;
t45 = pkin(4) * t34 + qJ(5) * t32;
t22 = -pkin(3) - t45;
t86 = t22 * t33;
t9 = t44 * qJD(4) - t32 * qJD(5);
t87 = t35 * t9;
t93 = (t86 + t88) * qJD(3) - qJD(4) * t8 - t87;
t89 = pkin(7) * t33;
t50 = pkin(3) * t35 + t89;
t20 = t50 * qJD(3) + qJD(2);
t39 = -t42 * qJD(4) + t34 * t20;
t92 = t45 * qJD(4) - t34 * qJD(5);
t91 = 0.2e1 * qJD(2);
t90 = 0.2e1 * qJD(5);
t85 = t32 * t36;
t84 = t35 * t36;
t82 = t28 + t30;
t29 = t33 ^ 2;
t31 = t35 ^ 2;
t81 = t29 - t31;
t80 = t29 + t31;
t79 = qJD(3) * t8;
t78 = qJD(3) * t34;
t77 = qJD(4) * t32;
t27 = qJD(4) * t34;
t76 = qJD(4) * t35;
t75 = qJD(4) * t36;
t74 = t33 * qJD(3);
t72 = t35 * qJD(3);
t71 = qJ(2) * qJD(3);
t70 = qJ(5) * qJD(3);
t69 = -0.2e1 * pkin(3) * qJD(4);
t59 = t36 * t72;
t68 = t32 * t20 + t21 * t27 + t34 * t59;
t67 = pkin(7) * t77;
t66 = pkin(7) * t27;
t65 = t32 * t76;
t64 = t32 * t75;
t63 = t34 * t76;
t62 = t32 * t74;
t61 = t32 * t27;
t60 = t33 * t72;
t58 = t35 * t70;
t57 = -pkin(4) + t85;
t56 = t82 * t35;
t54 = t81 * qJD(3);
t53 = 0.2e1 * t60;
t52 = t34 * t62;
t51 = t32 * t59;
t6 = qJ(5) * t33 + t42;
t7 = -t34 * t21 + t57 * t33;
t48 = t32 * t7 + t34 * t6;
t47 = t32 * t6 - t34 * t7;
t5 = t92 * t35 - t41 * t74;
t40 = -t5 + (t22 * t35 - t89) * qJD(4);
t1 = t58 + (qJD(5) - t64) * t33 + t68;
t2 = t57 * t72 - t39;
t37 = -t47 * qJD(4) + t1 * t34 + t2 * t32;
t17 = -t62 + t63;
t16 = t33 * t27 + t32 * t72;
t15 = t80 * t27;
t14 = -t34 * t74 - t65;
t13 = t33 * t77 - t34 * t72;
t12 = t80 * t77;
t4 = t39 - t51;
t3 = t33 * t64 - t68;
t10 = [0, 0, 0, 0, t91, qJ(2) * t91, -0.2e1 * t60, 0.2e1 * t54, 0, 0, 0, 0.2e1 * qJD(2) * t33 + 0.2e1 * t35 * t71, 0.2e1 * qJD(2) * t35 - 0.2e1 * t33 * t71, -0.2e1 * t30 * t60 - 0.2e1 * t31 * t61, 0.2e1 * t31 * t55 + 0.4e1 * t35 * t52, -0.2e1 * t33 * t65 - 0.2e1 * t81 * t78, 0.2e1 * t32 * t54 - 0.2e1 * t33 * t63, t53, 0.2e1 * (t21 * t72 - t31 * t75) * t34 + 0.2e1 * (t4 + t51) * t33, 0.2e1 * t31 * t64 + 0.2e1 * t3 * t33 + 0.2e1 * (-t19 + t24) * t72, 0.2e1 * (-t32 * t79 - t2) * t33 + 0.2e1 * (-qJD(3) * t7 + t8 * t27 + t32 * t5) * t35, 0.2e1 * t47 * t74 + 0.2e1 * (-t48 * qJD(4) - t1 * t32 + t2 * t34) * t35, 0.2e1 * (t8 * t78 + t1) * t33 + 0.2e1 * (qJD(3) * t6 - t34 * t5 + t8 * t77) * t35, 0.2e1 * t1 * t6 + 0.2e1 * t2 * t7 + 0.2e1 * t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t12, -t15, 0, -t12, (t48 * qJD(3) - t5) * t35 + (t37 + t79) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t82) * t53; 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t72, 0, -t36 * t74, -t59, -t35 * t55 - t52, -0.4e1 * t35 * t61 + t83 * t74, t16, -t13, 0, (-t32 * t84 - t50 * t34) * qJD(4) + (t49 * t32 - t24) * qJD(3), (t50 * t32 - t34 * t84) * qJD(4) + (-t34 * t88 + (pkin(3) * t34 + t85) * t33) * qJD(3), -t93 * t32 + t40 * t34, t37, t40 * t32 + t93 * t34, t37 * pkin(7) + t22 * t5 + t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t72, 0, 0, 0, 0, 0, t14, -t17, t14, qJD(3) * t56, t17, -t87 + (pkin(7) * t56 + t86) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t61, -0.2e1 * t55, 0, 0, 0, t32 * t69, t34 * t69, 0.2e1 * t22 * t77 - 0.2e1 * t34 * t9, 0, -0.2e1 * t22 * t27 - 0.2e1 * t32 * t9, 0.2e1 * t22 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t17, t72, t4, t3, (0.2e1 * pkin(4) - t85) * t72 + t39, (pkin(4) * t74 - qJ(5) * t76) * t34 + (t33 * t70 + (pkin(4) * qJD(4) - qJD(5)) * t35) * t32, 0.2e1 * t58 + (t90 - t64) * t33 + t68, -pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t13, -t16, 0, -t13, -t33 * t92 - t44 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t77, 0, -t66, t67, -t66, -t92, -t67, -t92 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, qJ(5) * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t14, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
