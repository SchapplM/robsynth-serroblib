% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRR14_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:37
% EndTime: 2024-09-27 18:43:39
% DurationCPUTime: 1.24s
% Computational Cost: add. (2435->167), mult. (7039->254), div. (0->0), fcn. (6307->10), ass. (0->153)
t120 = sin(qJ(4));
t121 = sin(qJ(3));
t124 = cos(qJ(4));
t125 = cos(qJ(3));
t197 = -t120 * t121 + t124 * t125;
t117 = sin(pkin(5));
t196 = t117 * (qJD(3) + qJD(4));
t118 = cos(pkin(5));
t195 = -0.2e1 * t118;
t194 = 0.2e1 * t118;
t65 = t197 * t196;
t193 = pkin(10) * t65;
t90 = t197 * t117;
t192 = t90 * pkin(4);
t191 = pkin(3) * t125;
t115 = t118 * pkin(3);
t114 = t118 * pkin(4);
t119 = sin(qJ(5));
t123 = cos(qJ(5));
t132 = t120 * t125 + t121 * t124;
t91 = t132 * t117;
t62 = t119 * t91 - t123 * t90;
t66 = t132 * t196;
t28 = -t62 * qJD(5) - t119 * t66 + t123 * t65;
t122 = sin(qJ(2));
t186 = pkin(1) * qJD(2);
t155 = t122 * t186;
t106 = t117 * t155;
t172 = qJD(3) * t121;
t145 = t117 * t172;
t105 = pkin(3) * t145;
t55 = pkin(4) * t66 + t105;
t45 = t106 + t55;
t63 = t119 * t90 + t123 * t91;
t126 = cos(qJ(2));
t113 = pkin(1) * t126 + pkin(2);
t95 = (-t113 - t191) * t117;
t69 = t95 - t192;
t190 = t69 * t28 + t45 * t63;
t97 = (-pkin(2) - t191) * t117;
t76 = t97 - t192;
t189 = t76 * t28 + t55 * t63;
t92 = t106 + t105;
t188 = t95 * t65 + t92 * t91;
t187 = t91 * t105 + t97 * t65;
t112 = pkin(3) * t124 + pkin(4);
t169 = qJD(4) * t124;
t152 = pkin(3) * t169;
t167 = qJD(5) * t123;
t73 = -t112 * t167 - t123 * t152 + (qJD(4) + qJD(5)) * t120 * t119 * pkin(3);
t185 = t118 * t73;
t177 = t117 * t125;
t111 = pkin(9) * t177;
t103 = pkin(1) * t122 + pkin(8) * t117;
t176 = t118 * t121;
t148 = t113 * t176;
t130 = -t103 * t125 - t148;
t75 = t111 - t130;
t184 = t120 * t75;
t156 = pkin(2) * t176;
t129 = -pkin(8) * t177 - t156;
t86 = t111 - t129;
t183 = t120 * t86;
t180 = t124 * t75;
t143 = -pkin(9) * t117 - t103;
t175 = t118 * t125;
t70 = t113 * t175 + t143 * t121 + t115;
t134 = -t120 * t70 - t180;
t88 = t90 * pkin(10);
t40 = -t134 + t88;
t182 = t123 * t40;
t179 = t124 * t86;
t149 = t117 * (-pkin(8) - pkin(9));
t135 = t121 * t149;
t83 = pkin(2) * t175 + t115 + t135;
t133 = -t120 * t83 - t179;
t43 = -t133 + t88;
t181 = t123 * t43;
t171 = qJD(3) * t125;
t146 = t118 * t171;
t154 = t126 * t186;
t178 = -t113 * t146 - t125 * t154;
t170 = qJD(4) * t120;
t168 = qJD(5) * t119;
t137 = t118 * t155;
t51 = (t143 * qJD(3) - t137) * t121 - t178;
t128 = (-t121 * t126 - t122 * t175) * t186;
t52 = t128 + (t143 * t125 - t148) * qJD(3);
t160 = -t120 * t52 - t124 * t51 - t70 * t169;
t14 = t75 * t170 + t160;
t64 = t66 * pkin(10);
t11 = -t14 - t64;
t140 = -t120 * t51 + t124 * t52;
t15 = t134 * qJD(4) + t140;
t12 = t15 - t193;
t144 = -pkin(10) * t91 + t114;
t37 = t124 * t70 + t144 - t184;
t166 = -t123 * t11 - t119 * t12 - t37 * t167;
t29 = t63 * qJD(5) + t119 * t65 + t123 * t66;
t142 = -t119 * t11 + t123 * t12;
t3 = (-t119 * t37 - t182) * qJD(5) + t142;
t165 = t3 * t118 + t69 * t29 + t45 * t62;
t107 = pkin(2) * t146;
t84 = qJD(3) * t135 + t107;
t85 = (t125 * t149 - t156) * qJD(3);
t159 = -t120 * t85 - t124 * t84 - t83 * t169;
t32 = t86 * t170 + t159;
t24 = -t32 - t64;
t139 = -t120 * t84 + t124 * t85;
t33 = t133 * qJD(4) + t139;
t25 = t33 - t193;
t141 = -t119 * t24 + t123 * t25;
t42 = t124 * t83 + t144 - t183;
t6 = (-t119 * t42 - t181) * qJD(5) + t141;
t164 = t6 * t118 + t76 * t29 + t55 * t62;
t163 = t15 * t118 + t95 * t66 - t92 * t90;
t162 = -t119 * t25 - t123 * t24 - t42 * t167;
t161 = -t90 * t105 + t33 * t118 + t97 * t66;
t158 = t124 * t115;
t157 = t123 * t114;
t153 = pkin(3) * t170;
t151 = pkin(4) * t168;
t150 = pkin(4) * t167;
t116 = t117 ^ 2;
t147 = t116 * t171;
t110 = t117 * t171;
t138 = qJD(3) * (-pkin(2) - t113);
t136 = t125 * t155;
t2 = t40 * t168 + t166;
t5 = t43 * t168 + t162;
t127 = (-t120 * t167 + (-t119 * t124 - t120 * t123) * qJD(4)) * pkin(3);
t101 = 0.2e1 * t121 * t147;
t100 = t116 * t121 * t155;
t99 = t110 * t194;
t98 = t145 * t195;
t94 = t129 * qJD(3);
t93 = pkin(8) * t145 - t107;
t89 = 0.2e1 * (-t121 ^ 2 + t125 ^ 2) * t116 * qJD(3);
t87 = t94 * t118;
t74 = -t112 * t168 + t127;
t68 = t74 * t118;
t60 = t66 * t195;
t59 = t65 * t194;
t58 = t130 * qJD(3) + t128;
t57 = (qJD(3) * t103 + t137) * t121 + t178;
t56 = t58 * t118;
t44 = 0.2e1 * t91 * t65;
t30 = 0.2e1 * t65 * t90 - 0.2e1 * t66 * t91;
t27 = t29 * t195;
t26 = t28 * t194;
t16 = 0.2e1 * t63 * t28;
t7 = -0.2e1 * t28 * t62 - 0.2e1 * t29 * t63;
t1 = [0, 0, 0, 0, -0.2e1 * t155, -0.2e1 * t154, t101, t89, t99, t98, 0, 0.2e1 * t56 + 0.2e1 * (-t113 * t172 - t136) * t116, -0.2e1 * t113 * t147 + 0.2e1 * t118 * t57 + 0.2e1 * t100, t44, t30, t59, t60, 0, 0.2e1 * t163, 0.2e1 * t118 * t14 + 0.2e1 * t188, t16, t7, t26, t27, 0, 0.2e1 * t165, 0.2e1 * t118 * t2 + 0.2e1 * t190; 0, 0, 0, 0, -t155, -t154, t101, t89, t99, t98, 0, t56 + t87 + (t121 * t138 - t136) * t116, t100 + (t57 + t93) * t118 + t125 * t116 * t138, t44, t30, t59, t60, 0, t161 + t163, (t14 + t32) * t118 + t187 + t188, t16, t7, t26, t27, 0, t164 + t165, (t2 + t5) * t118 + t189 + t190; 0, 0, 0, 0, 0, 0, t101, t89, t99, t98, 0, -0.2e1 * pkin(2) * t116 * t172 + 0.2e1 * t87, -0.2e1 * pkin(2) * t147 + 0.2e1 * t118 * t93, t44, t30, t59, t60, 0, 0.2e1 * t161, 0.2e1 * t118 * t32 + 0.2e1 * t187, t16, t7, t26, t27, 0, 0.2e1 * t164, 0.2e1 * t118 * t5 + 0.2e1 * t189; 0, 0, 0, 0, 0, 0, 0, 0, t110, -t145, 0, t58, t57, 0, 0, t65, -t66, 0, (-t180 + (-t70 - t115) * t120) * qJD(4) + t140, (-t158 + t184) * qJD(4) + t160, 0, 0, t28, -t29, 0, t68 + t3, t2 + t185; 0, 0, 0, 0, 0, 0, 0, 0, t110, -t145, 0, t94, t93, 0, 0, t65, -t66, 0, (-t179 + (-t83 - t115) * t120) * qJD(4) + t139, (-t158 + t183) * qJD(4) + t159, 0, 0, t28, -t29, 0, t68 + t6, t5 + t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t153, -0.2e1 * t152, 0, 0, 0, 0, 0, 0.2e1 * t74, 0.2e1 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t66, 0, t15, t14, 0, 0, t28, -t29, 0, (-t182 + (-t37 - t114) * t119) * qJD(5) + t142, (t119 * t40 - t157) * qJD(5) + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t66, 0, t33, t32, 0, 0, t28, -t29, 0, (-t181 + (-t42 - t114) * t119) * qJD(5) + t141, (t119 * t43 - t157) * qJD(5) + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, -t152, 0, 0, 0, 0, 0, (-pkin(4) - t112) * t168 + t127, t73 - t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t151, -0.2e1 * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, 0, t3, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, -t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
