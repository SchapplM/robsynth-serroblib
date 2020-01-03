% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPRP3
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPRP3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:51:13
% EndTime: 2019-12-31 19:51:16
% DurationCPUTime: 0.74s
% Computational Cost: add. (848->96), mult. (1951->143), div. (0->0), fcn. (1673->6), ass. (0->73)
t73 = sin(pkin(8));
t74 = cos(pkin(8));
t102 = t73 ^ 2 + t74 ^ 2;
t131 = t102 * qJD(3);
t133 = 0.2e1 * t131;
t101 = pkin(1) * qJD(2);
t77 = cos(qJ(2));
t98 = t77 * t101;
t61 = qJD(3) + t98;
t132 = t102 * t61;
t76 = sin(qJ(2));
t65 = t76 * pkin(1) + qJ(3);
t70 = t74 * pkin(7);
t51 = t74 * t65 + t70;
t75 = sin(qJ(4));
t122 = -pkin(7) - t65;
t118 = cos(qJ(4));
t96 = t118 * t73;
t80 = t122 * t96;
t95 = t118 * t74;
t10 = -qJD(4) * t80 - t61 * t95 + (qJD(4) * t51 + t61 * t73) * t75;
t111 = t75 * t74;
t88 = qJD(4) * t118;
t97 = t75 * t122;
t11 = t51 * t88 + t61 * t111 + (qJD(4) * t97 + t118 * t61) * t73;
t35 = t75 * t51 - t80;
t36 = t118 * t51 + t73 * t97;
t100 = qJD(4) * t75;
t49 = t73 * t100 - t74 * t88;
t56 = t96 + t111;
t50 = t56 * qJD(4);
t55 = t75 * t73 - t95;
t125 = t10 * t55 + t11 * t56 - t35 * t49 - t36 * t50;
t130 = 0.2e1 * t125;
t60 = t74 * qJ(3) + t70;
t110 = -pkin(7) - qJ(3);
t92 = t110 * t73;
t79 = t118 * t92;
t87 = t118 * qJD(3);
t29 = -qJD(4) * t79 - t74 * t87 + (qJD(3) * t73 + qJD(4) * t60) * t75;
t30 = t60 * t88 + qJD(3) * t111 + (t110 * t100 + t87) * t73;
t41 = t75 * t60 - t79;
t42 = t118 * t60 + t75 * t92;
t126 = t29 * t55 + t30 * t56 - t41 * t49 - t42 * t50;
t129 = 0.2e1 * t126;
t128 = t126 + t125;
t127 = 0.2e1 * t49 * t55 - 0.2e1 * t56 * t50;
t124 = 2 * qJD(5);
t123 = t77 * pkin(1);
t27 = t50 * pkin(4) + t49 * qJ(5) - t56 * qJD(5);
t67 = t76 * t101;
t19 = t27 + t67;
t66 = -t74 * pkin(3) - pkin(2);
t37 = t55 * pkin(4) - t56 * qJ(5) + t66;
t34 = t37 - t123;
t120 = t19 * t55 + t34 * t50;
t119 = -t19 * t56 + t34 * t49;
t113 = t66 * t49;
t112 = t66 * t50;
t109 = t27 * t55 + t37 * t50;
t108 = -t27 * t56 + t37 * t49;
t59 = t66 - t123;
t106 = t59 * t50 + t55 * t67;
t105 = -t59 * t49 + t56 * t67;
t99 = 0.2e1 * t55 * t50;
t94 = -t36 * t10 + t35 * t11;
t93 = -t42 * t29 + t41 * t30;
t83 = t73 * t67;
t82 = t74 * t67;
t78 = -t10 * t42 + t11 * t41 - t36 * t29 + t35 * t30;
t40 = -0.2e1 * t56 * t49;
t28 = pkin(4) * t49 - t50 * qJ(5) - t55 * qJD(5);
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t67, -0.2e1 * t98, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t82, 0.2e1 * t83, 0.2e1 * t132, 0.2e1 * (-pkin(2) - t123) * t67 + 0.2e1 * t65 * t132, t40, t127, 0, t99, 0, 0, 0.2e1 * t106, 0.2e1 * t105, t130, 0.2e1 * t59 * t67 + 0.2e1 * t94, t40, 0, -t127, 0, 0, t99, 0.2e1 * t120, t130, 0.2e1 * t119, 0.2e1 * t34 * t19 + 0.2e1 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t98, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t83, t131 + t132, -pkin(2) * t67 + qJ(3) * t132 + t131 * t65, t40, t127, 0, t99, 0, 0, t106 + t112, t105 - t113, t128, t66 * t67 + t78, t40, 0, -t127, 0, 0, t99, t109 + t120, t128, t108 + t119, t19 * t37 + t34 * t27 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, qJ(3) * t133, t40, t127, 0, t99, 0, 0, 0.2e1 * t112, -0.2e1 * t113, t129, 0.2e1 * t93, t40, 0, -t127, 0, 0, t99, 0.2e1 * t109, t129, 0.2e1 * t108, 0.2e1 * t37 * t27 + 0.2e1 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, t50, -t49, 0, t67, 0, 0, 0, 0, 0, 0, t50, 0, t49, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t49, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t49, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, -t50, 0, -t11, t10, 0, 0, 0, -t49, 0, 0, t50, 0, -t11, t28, -t10, -t11 * pkin(4) - t10 * qJ(5) + t36 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, -t50, 0, -t30, t29, 0, 0, 0, -t49, 0, 0, t50, 0, -t30, t28, -t29, -t30 * pkin(4) - t29 * qJ(5) + t42 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, qJ(5) * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
