% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:40
% EndTime: 2019-12-31 19:01:44
% DurationCPUTime: 1.18s
% Computational Cost: add. (1764->142), mult. (3774->261), div. (0->0), fcn. (3346->8), ass. (0->101)
t60 = sin(qJ(5));
t58 = t60 ^ 2;
t63 = cos(qJ(5));
t59 = t63 ^ 2;
t124 = (t58 - t59) * qJD(5);
t61 = sin(qJ(4));
t62 = sin(qJ(3));
t107 = t61 * t62;
t123 = qJD(3) + qJD(4);
t64 = cos(qJ(3));
t116 = cos(qJ(4));
t81 = qJD(4) * t116;
t84 = t116 * t64;
t32 = -qJD(3) * t84 + t123 * t107 - t64 * t81;
t106 = t61 * t64;
t43 = t116 * t62 + t106;
t33 = t123 * t43;
t42 = -t84 + t107;
t127 = t42 * t32 - t43 * t33;
t126 = t58 + t59;
t52 = sin(pkin(9)) * pkin(1) + pkin(6);
t118 = pkin(7) + t52;
t39 = t118 * t64;
t28 = -t118 * t107 + t116 * t39;
t53 = -cos(pkin(9)) * pkin(1) - pkin(2);
t46 = -t64 * pkin(3) + t53;
t67 = -t42 * pkin(4) + t43 * pkin(8) - t46;
t66 = t63 * t67;
t7 = -t60 * t28 - t66;
t8 = t63 * t28 - t60 * t67;
t125 = -t60 * t7 + t63 * t8;
t122 = 0.2e1 * qJD(3);
t121 = t33 * pkin(4);
t78 = t118 * t116;
t27 = t61 * t39 + t62 * t78;
t57 = qJD(5) * t63;
t82 = qJD(3) * t118;
t38 = t62 * t82;
t95 = t64 * qJD(3);
t10 = t28 * qJD(4) - t61 * t38 + t78 * t95;
t6 = t10 * t60;
t117 = t27 * t57 + t6;
t115 = t10 * t43;
t114 = t27 * t10;
t113 = t27 * t61;
t112 = t42 * t33;
t111 = t42 * t61;
t110 = t43 * t32;
t109 = t58 * t32;
t30 = t59 * t32;
t108 = t60 * t33;
t105 = t63 * t32;
t104 = t63 * t33;
t103 = t43 * t104 - t42 * t105;
t56 = -t116 * pkin(3) - pkin(4);
t99 = pkin(3) * qJD(4);
t88 = t61 * t99;
t101 = t56 * t57 + t60 * t88;
t97 = qJD(5) * t60;
t96 = t62 * qJD(3);
t94 = 0.2e1 * t112;
t93 = t60 * t105;
t92 = t53 * t122;
t91 = pkin(4) * t97;
t90 = pkin(4) * t57;
t89 = pkin(3) * t96;
t87 = t43 * t97;
t86 = t60 * t57;
t85 = t62 * t95;
t83 = t126 * t32;
t40 = t43 ^ 2;
t80 = t40 * t86;
t79 = pkin(3) * t81;
t77 = t60 * t8 + t63 * t7;
t75 = t10 * t42 + t27 * t33;
t55 = t61 * pkin(3) + pkin(8);
t73 = t42 * t55 - t43 * t56;
t72 = t56 * t97 - t63 * t88;
t71 = -t60 * t32 + t43 * t57;
t15 = t87 + t105;
t14 = t42 * t97 - t104;
t70 = t126 * t116;
t69 = (-t116 * t42 + t43 * t61) * qJD(4);
t68 = t32 * pkin(8) + t121 + t89;
t9 = t27 * qJD(4) + t82 * t106 + t116 * t38;
t2 = qJD(5) * t66 + t28 * t97 - t60 * t68 + t63 * t9;
t3 = -t8 * qJD(5) + t60 * t9 + t63 * t68;
t1 = -t77 * qJD(5) - t2 * t63 - t3 * t60;
t65 = pkin(3) * t69 - t32 * t56 - t33 * t55;
t49 = -0.2e1 * t86;
t48 = 0.2e1 * t86;
t41 = -0.2e1 * t124;
t37 = t70 * t99;
t25 = t27 * t97;
t19 = t43 * t30;
t18 = t43 * t109;
t16 = t42 * t57 + t108;
t13 = -t30 - t109;
t12 = t43 * t124 + t93;
t4 = -0.4e1 * t43 * t86 + t109 - t30;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t85, (-t62 ^ 2 + t64 ^ 2) * t122, 0, -0.2e1 * t85, 0, 0, t62 * t92, t64 * t92, 0, 0, -0.2e1 * t110, 0.2e1 * t127, 0, t94, 0, 0, 0.2e1 * t46 * t33 + 0.2e1 * t42 * t89, -0.2e1 * t46 * t32 + 0.2e1 * t43 * t89, -0.2e1 * t27 * t32 - 0.2e1 * t28 * t33 + 0.2e1 * t9 * t42 + 0.2e1 * t115, -0.2e1 * t28 * t9 + 0.2e1 * t46 * t89 + 0.2e1 * t114, -0.2e1 * t19 - 0.2e1 * t80, 0.2e1 * t40 * t124 + 0.4e1 * t43 * t93, -0.2e1 * t42 * t87 + 0.2e1 * t103, -0.2e1 * t18 + 0.2e1 * t80, -0.2e1 * t43 * t108 - 0.2e1 * t71 * t42, t94, 0.2e1 * t71 * t27 + 0.2e1 * t3 * t42 + 0.2e1 * t7 * t33 + 0.2e1 * t43 * t6, 0.2e1 * t63 * t115 - 0.2e1 * t15 * t27 + 0.2e1 * t2 * t42 - 0.2e1 * t8 * t33, 0.2e1 * t77 * t32 + 0.2e1 * (-qJD(5) * t125 + t2 * t60 - t3 * t63) * t43, -0.2e1 * t8 * t2 + 0.2e1 * t7 * t3 + 0.2e1 * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28 * t32 - t9 * t43 + t75, 0, 0, 0, 0, 0, 0, 0, t127 * t63 + t103, 0, t1 * t43 - t125 * t32 + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t110 + 0.2e1 * t112, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t18 - 0.2e1 * t19 + 0.2e1 * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, -t96, 0, -t52 * t95, t52 * t96, 0, 0, 0, 0, -t32, 0, -t33, 0, -t10, t9, (t116 * t32 - t33 * t61 + t69) * pkin(3), (-t116 * t10 - t61 * t9 + (t116 * t28 + t113) * qJD(4)) * pkin(3), -t12, t4, t16, t12, -t14, 0, t25 + (-t73 * qJD(5) - t10) * t63 + t65 * t60, t63 * t65 + t73 * t97 + t117, t1, t10 * t56 + (t125 * t116 + t113) * t99 + t1 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t95, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, 0, (-t116 * t33 - t32 * t61 + (t116 * t43 + t111) * qJD(4)) * pkin(3), 0, 0, 0, 0, 0, 0, t14, t16, t13, t33 * t56 - t55 * t83 + (t43 * t70 + t111) * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t88, -0.2e1 * t79, 0, 0, t48, t41, 0, t49, 0, 0, 0.2e1 * t72, 0.2e1 * t101, 0.2e1 * t37, 0.2e1 * (t55 * t70 + t56 * t61) * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, 0, -t33, 0, -t10, t9, 0, 0, -t12, t4, t16, t12, -t14, 0, t25 + (pkin(4) * t32 - pkin(8) * t33) * t60 + (-t10 + (-pkin(4) * t43 - pkin(8) * t42) * qJD(5)) * t63, pkin(4) * t15 + pkin(8) * t14 + t117, t1, -t10 * pkin(4) + pkin(8) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t32, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, t13, -pkin(8) * t83 - t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t79, 0, 0, t48, t41, 0, t49, 0, 0, t72 - t91, -t90 + t101, t37, (-pkin(4) * t61 + pkin(8) * t70) * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t41, 0, t49, 0, 0, -0.2e1 * t91, -0.2e1 * t90, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, -t71, t33, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, -t97, 0, -t55 * t57 - t60 * t79, t55 * t97 - t63 * t79, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, -t97, 0, -pkin(8) * t57, pkin(8) * t97, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
