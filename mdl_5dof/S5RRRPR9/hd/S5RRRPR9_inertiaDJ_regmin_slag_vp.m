% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:23:55
% EndTime: 2021-01-15 23:24:05
% DurationCPUTime: 1.56s
% Computational Cost: add. (1881->202), mult. (4851->391), div. (0->0), fcn. (4233->8), ass. (0->116)
t104 = sin(qJ(2));
t103 = sin(qJ(3));
t107 = cos(qJ(2));
t135 = t107 * qJD(2);
t125 = t103 * t135;
t106 = cos(qJ(3));
t138 = qJD(3) * t106;
t152 = t104 * t138 + t125;
t151 = -0.4e1 * t104;
t100 = sin(pkin(9));
t101 = cos(pkin(9));
t139 = qJD(3) * t103;
t61 = -t100 * t139 + t101 * t138;
t97 = t104 ^ 2;
t119 = (-t107 ^ 2 + t97) * qJD(2);
t98 = t106 ^ 2;
t144 = t103 ^ 2 - t98;
t120 = t144 * qJD(3);
t150 = pkin(3) * t100;
t149 = pkin(6) * t103;
t148 = qJ(4) + pkin(7);
t136 = t106 * qJD(4);
t140 = t106 * t107;
t94 = t104 * qJD(2);
t126 = t103 * t94;
t115 = pkin(2) * t104 - pkin(7) * t107;
t72 = t115 * qJD(2);
t146 = pkin(6) * t126 + t106 * t72;
t116 = -t107 * pkin(2) - t104 * pkin(7);
t78 = -pkin(1) + t116;
t88 = pkin(6) * t140;
t16 = -t104 * t136 + (pkin(3) * t104 - qJ(4) * t140) * qJD(2) + (-t88 + (qJ(4) * t104 - t78) * t103) * qJD(3) + t146;
t141 = t104 * t106;
t147 = -t103 * t72 - t78 * t138;
t20 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t141 + (-qJD(4) * t104 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t107) * t103 - t147;
t9 = t100 * t16 + t101 * t20;
t68 = t106 * t78;
t41 = -qJ(4) * t141 + t68 + (-pkin(3) - t149) * t107;
t142 = t103 * t104;
t145 = t103 * t78 + t88;
t44 = -qJ(4) * t142 + t145;
t22 = t100 * t41 + t101 * t44;
t79 = t148 * t103;
t80 = t148 * t106;
t46 = -t100 * t79 + t101 * t80;
t73 = pkin(3) * t142 + t104 * pkin(6);
t137 = qJD(3) * t107;
t134 = -0.2e1 * pkin(1) * qJD(2);
t133 = -0.2e1 * pkin(2) * qJD(3);
t92 = pkin(6) * t135;
t48 = t152 * pkin(3) + t92;
t132 = pkin(3) * t94;
t93 = pkin(3) * t139;
t91 = -t106 * pkin(3) - pkin(2);
t131 = t103 * t137;
t129 = t106 * t137;
t124 = t103 * t138;
t123 = t104 * t135;
t122 = t106 * t135;
t8 = -t100 * t20 + t101 * t16;
t21 = -t100 * t44 + t101 * t41;
t121 = qJD(3) * t148;
t57 = -t103 * t121 + t136;
t58 = -t103 * qJD(4) - t106 * t121;
t33 = -t100 * t57 + t101 * t58;
t45 = -t100 * t80 - t101 * t79;
t118 = 0.2e1 * t123;
t117 = t103 * t122;
t34 = t100 * t58 + t101 * t57;
t102 = sin(qJ(5));
t105 = cos(qJ(5));
t56 = -t100 * t142 + t101 * t141;
t12 = -t107 * pkin(4) - t56 * pkin(8) + t21;
t66 = t100 * t106 + t101 * t103;
t55 = t66 * t104;
t13 = -t55 * pkin(8) + t22;
t114 = t102 * t13 - t105 * t12;
t113 = t102 * t12 + t105 * t13;
t35 = -t66 * pkin(8) + t45;
t65 = t100 * t103 - t101 * t106;
t36 = -t65 * pkin(8) + t46;
t112 = t102 * t36 - t105 * t35;
t111 = t102 * t35 + t105 * t36;
t29 = t102 * t56 + t105 * t55;
t30 = -t102 * t55 + t105 * t56;
t39 = t102 * t66 + t105 * t65;
t40 = -t102 * t65 + t105 * t66;
t90 = t101 * pkin(3) + pkin(4);
t110 = t102 * t90 + t105 * t150;
t109 = t102 * t150 - t105 * t90;
t60 = t66 * qJD(3);
t108 = t106 * t94 + t131;
t32 = t100 * t125 - t101 * t122 + t104 * t60;
t6 = pkin(4) * t94 + t32 * pkin(8) + t8;
t31 = -t100 * t122 - t101 * t125 - t61 * t104;
t7 = t31 * pkin(8) + t9;
t2 = -t113 * qJD(5) - t102 * t7 + t105 * t6;
t1 = t114 * qJD(5) - t102 * t6 - t105 * t7;
t85 = -0.2e1 * t123;
t53 = t110 * qJD(5);
t52 = t109 * qJD(5);
t50 = t65 * pkin(4) + t91;
t47 = t60 * pkin(4) + t93;
t42 = t55 * pkin(4) + t73;
t27 = -t145 * qJD(3) + t146;
t26 = t108 * pkin(6) + t147;
t25 = -t60 * pkin(8) + t34;
t24 = -t61 * pkin(8) + t33;
t23 = -t31 * pkin(4) + t48;
t18 = t40 * qJD(5) + t102 * t61 + t105 * t60;
t17 = -t39 * qJD(5) - t102 * t60 + t105 * t61;
t11 = t30 * qJD(5) - t102 * t32 - t105 * t31;
t10 = -t29 * qJD(5) + t102 * t31 - t105 * t32;
t4 = -t111 * qJD(5) - t102 * t25 + t105 * t24;
t3 = t112 * qJD(5) - t102 * t24 - t105 * t25;
t5 = [0, 0, 0, t118, -0.2e1 * t119, 0, 0, 0, t104 * t134, t107 * t134, 0.2e1 * t98 * t123 - 0.2e1 * t97 * t124, t117 * t151 + 0.2e1 * t97 * t120, 0.2e1 * t104 * t131 + 0.2e1 * t106 * t119, -0.2e1 * t103 * t119 + 0.2e1 * t104 * t129, t85, 0.2e1 * t68 * t94 - 0.2e1 * t27 * t107 + 0.2e1 * (t103 * t123 + t97 * t138) * pkin(6), -0.2e1 * t26 * t107 - 0.2e1 * t145 * t94 + 0.2e1 * (t106 * t118 - t97 * t139) * pkin(6), -0.2e1 * t8 * t107 + 0.2e1 * t21 * t94 - 0.2e1 * t73 * t31 + 0.2e1 * t48 * t55, 0.2e1 * t9 * t107 - 0.2e1 * t22 * t94 - 0.2e1 * t73 * t32 + 0.2e1 * t48 * t56, 0.2e1 * t21 * t32 + 0.2e1 * t22 * t31 - 0.2e1 * t9 * t55 - 0.2e1 * t8 * t56, 0.2e1 * t21 * t8 + 0.2e1 * t22 * t9 + 0.2e1 * t73 * t48, 0.2e1 * t30 * t10, -0.2e1 * t10 * t29 - 0.2e1 * t30 * t11, -0.2e1 * t10 * t107 + 0.2e1 * t30 * t94, 0.2e1 * t11 * t107 - 0.2e1 * t29 * t94, t85, -0.2e1 * t2 * t107 + 0.2e1 * t42 * t11 - 0.2e1 * t114 * t94 + 0.2e1 * t23 * t29, -0.2e1 * t1 * t107 + 0.2e1 * t42 * t10 - 0.2e1 * t113 * t94 + 0.2e1 * t23 * t30; 0, 0, 0, 0, 0, t135, -t94, 0, -t92, pkin(6) * t94, -t104 * t120 + t117, t124 * t151 - t144 * t135, t126 - t129, t108, 0, (pkin(7) * t140 + (-pkin(2) * t106 + t149) * t104) * qJD(3) + (t116 * t103 - t88) * qJD(2), (pkin(6) * t141 + t115 * t103) * qJD(3) + (t116 * t106 + t107 * t149) * qJD(2), -t33 * t107 - t91 * t31 + t45 * t94 + t48 * t65 + t55 * t93 + t73 * t60, t34 * t107 - t91 * t32 - t46 * t94 + t48 * t66 + t56 * t93 + t73 * t61, -t21 * t61 - t22 * t60 + t46 * t31 + t45 * t32 - t33 * t56 - t34 * t55 - t9 * t65 - t8 * t66, t21 * t33 + t22 * t34 + t8 * t45 + t9 * t46 + t48 * t91 + t73 * t93, t10 * t40 + t30 * t17, -t10 * t39 - t40 * t11 - t17 * t29 - t30 * t18, -t17 * t107 + t40 * t94, t18 * t107 - t39 * t94, 0, -t4 * t107 + t50 * t11 - t112 * t94 + t42 * t18 + t23 * t39 + t47 * t29, t50 * t10 - t3 * t107 - t111 * t94 + t42 * t17 + t23 * t40 + t47 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t124, -0.2e1 * t120, 0, 0, 0, t103 * t133, t106 * t133, 0.2e1 * t91 * t60 + 0.2e1 * t65 * t93, 0.2e1 * t91 * t61 + 0.2e1 * t66 * t93, -0.2e1 * t33 * t66 - 0.2e1 * t34 * t65 - 0.2e1 * t45 * t61 - 0.2e1 * t46 * t60, 0.2e1 * t45 * t33 + 0.2e1 * t46 * t34 + 0.2e1 * t91 * t93, 0.2e1 * t40 * t17, -0.2e1 * t17 * t39 - 0.2e1 * t40 * t18, 0, 0, 0, 0.2e1 * t50 * t18 + 0.2e1 * t47 * t39, 0.2e1 * t50 * t17 + 0.2e1 * t47 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104 * t139 + t122, -t152, t94, t27, t26, t101 * t132 + t8, -t100 * t132 - t9, (t100 * t31 + t101 * t32) * pkin(3), (t100 * t9 + t101 * t8) * pkin(3), 0, 0, t10, -t11, t94, t53 * t107 - t109 * t94 + t2, -t52 * t107 - t110 * t94 + t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, -t139, 0, -pkin(7) * t138, pkin(7) * t139, t33, -t34, (-t100 * t60 - t101 * t61) * pkin(3), (t100 * t34 + t101 * t33) * pkin(3), 0, 0, t17, -t18, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t53, 0.2e1 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t32, 0, t48, 0, 0, 0, 0, 0, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t61, 0, t93, 0, 0, 0, 0, 0, t18, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, t94, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t18, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
