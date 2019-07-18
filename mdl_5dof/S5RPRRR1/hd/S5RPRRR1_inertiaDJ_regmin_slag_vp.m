% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_inertiaDJ_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:02
% EndTime: 2019-07-18 13:26:06
% DurationCPUTime: 0.81s
% Computational Cost: add. (272->113), mult. (1107->241), div. (0->0), fcn. (928->6), ass. (0->98)
t36 = sin(qJ(3));
t112 = -0.4e1 * t36;
t35 = sin(qJ(4));
t29 = t35 ^ 2;
t39 = cos(qJ(3));
t27 = t39 * qJD(3);
t37 = cos(qJ(5));
t34 = sin(qJ(5));
t88 = qJD(5) * t34;
t44 = t37 * t27 - t36 * t88;
t111 = t44 * t29;
t30 = t36 ^ 2;
t33 = t39 ^ 2;
t96 = t30 + t33;
t110 = qJ(2) * qJD(4) * t96;
t84 = qJ(2) * qJD(3);
t94 = qJD(2) * t36;
t109 = -t39 * t84 - t94;
t38 = cos(qJ(4));
t32 = t38 ^ 2;
t98 = t29 - t32;
t108 = qJD(4) * t98;
t31 = t37 ^ 2;
t99 = t34 ^ 2 - t31;
t61 = t99 * qJD(5);
t107 = -0.4e1 * t35;
t106 = 0.2e1 * qJD(2);
t100 = t39 * t37;
t92 = qJD(3) * t38;
t51 = -qJD(5) + t92;
t86 = qJD(5) * t38;
t52 = qJD(3) - t86;
t91 = qJD(4) * t35;
t80 = t37 * t91;
t4 = t51 * t100 + (t52 * t34 - t80) * t36;
t105 = t4 * t34;
t104 = t4 * t37;
t103 = t35 * t36;
t102 = t36 * t38;
t101 = t39 * t34;
t97 = t30 - t33;
t95 = qJ(2) * t39;
t93 = qJD(2) * t39;
t90 = qJD(4) * t38;
t89 = qJD(4) * t39;
t87 = qJD(5) * t37;
t85 = qJD(5) * t39;
t26 = t36 * qJD(3);
t82 = t29 * t93;
t81 = t36 * t91;
t79 = t35 * t89;
t78 = t38 * t89;
t77 = t36 * t87;
t76 = t35 * t87;
t75 = t34 * t87;
t74 = t35 * t90;
t73 = t35 * t93;
t72 = t36 * t27;
t71 = t37 * t90;
t70 = t38 * t27;
t69 = t34 * t86;
t68 = t37 * t86;
t67 = t38 * t93;
t66 = t34 * t85;
t65 = t37 * t85;
t63 = qJD(2) * t96;
t59 = t97 * qJD(3);
t58 = 0.2e1 * t74;
t57 = t34 * t81;
t56 = t36 * t80;
t55 = t34 * t71;
t54 = t35 * t70;
t53 = t30 * t74;
t50 = t39 * t58;
t16 = t34 * t102 + t100;
t17 = t37 * t102 - t101;
t49 = -t16 * t37 - t17 * t34;
t48 = -t38 * t101 + t36 * t37;
t42 = t34 * t27 + t77;
t43 = t37 * t26 + t66;
t3 = t42 * t38 - t43 - t57;
t47 = -t16 * t91 + t38 * t3;
t46 = t70 - t81;
t13 = t38 * t26 + t79;
t14 = t35 * t27 + t36 * t90;
t45 = t35 * t26 - t78;
t41 = -t35 * t88 + t71;
t40 = -t29 * t42 - 0.2e1 * t38 * t57;
t12 = t69 + t80;
t11 = t34 * t91 - t68;
t10 = -t34 * t90 - t76;
t8 = (t38 * t100 + t34 * t36) * qJ(2);
t7 = t48 * qJ(2);
t6 = t17 * t91;
t5 = t108 * t36 - t54;
t2 = t48 * qJD(2) + (t52 * t100 + (t51 * t36 + t79) * t34) * qJ(2);
t1 = -t37 * t67 + t109 * t34 + (t13 * t37 + t38 * t66 - t77) * qJ(2);
t9 = [0, 0, 0, 0, t106, qJ(2) * t106, 0.2e1 * t72, -0.2e1 * t59, 0, 0, 0, 0, 0, 0.2e1 * t32 * t72 - 0.2e1 * t53, 0.2e1 * t108 * t30 + t54 * t112, 0.2e1 * t36 * t79 + 0.2e1 * t97 * t92, -0.2e1 * t35 * t59 + 0.2e1 * t36 * t78, -0.2e1 * t72, 0.2e1 * t38 * t110 + 0.2e1 * t35 * t63, -0.2e1 * t35 * t110 + 0.2e1 * t38 * t63, 0.2e1 * t17 * t4, -0.2e1 * t4 * t16 - 0.2e1 * t17 * t3, 0.2e1 * t4 * t103 + 0.2e1 * t14 * t17, -0.2e1 * t3 * t103 - 0.2e1 * t14 * t16, 0.2e1 * t29 * t72 + 0.2e1 * t53, 0.2e1 * (t16 * t95 + t36 * t7) * t90 + 0.2e1 * ((-t16 * t84 + t2) * t36 + (qJ(2) * t3 + qJD(2) * t16 + qJD(3) * t7) * t39) * t35, 0.2e1 * (t17 * t95 - t36 * t8) * t90 + 0.2e1 * ((-t17 * t84 + t1) * t36 + (qJ(2) * t4 + qJD(2) * t17 - qJD(3) * t8) * t39) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t27, 0, 0, 0, 0, 0, t13, -t45, 0, 0, 0, 0, 0, t40 - t47, t6 + (-t4 - 0.2e1 * t56) * t38 - t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, 0, t109, t36 * t84 - t93, -t5, t74 * t112 - t98 * t27, t45, t13, 0, -t46 * qJ(2) - t38 * t94, t14 * qJ(2) + t35 * t94, t35 * t104 + t17 * t41, t49 * t90 + (-t3 * t37 - t105 + (t16 * t34 - t17 * t37) * qJD(5)) * t35, t6 + (-t4 + 0.2e1 * t56) * t38 + t111, t40 + t47, t5, t34 * t82 + t7 * t91 - t2 * t38 + (t34 * t50 + (-t34 * t26 + t65) * t29) * qJ(2), t37 * t82 - t8 * t91 - t1 * t38 + (-t29 * t43 + t37 * t50) * qJ(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -0.2e1 * t108, 0, 0, 0, 0, 0, -0.2e1 * t29 * t75 + 0.2e1 * t31 * t74, t55 * t107 + 0.2e1 * t29 * t61, 0.2e1 * t108 * t37 + 0.2e1 * t35 * t69, -0.2e1 * t108 * t34 + 0.2e1 * t35 * t68, -0.2e1 * t74, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t14, t26, t45 * qJ(2) - t73, t13 * qJ(2) - t67, t17 * t87 + t105, t49 * qJD(5) - t34 * t3 + t104, t14 * t34 + t36 * t76, t44 * t35 + t36 * t71, 0, -t37 * t73 + (t35 * t43 - t39 * t71) * qJ(2), t34 * t73 + (-t34 * t45 + t35 * t65) * qJ(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t90, 0, 0, 0, 0, 0, -t12, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t91, 0, 0, 0, -t35 * t61 + t55, t75 * t107 - t99 * t90, t11, t12, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t75, -0.2e1 * t61, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, t14, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t10, t91, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t88, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
