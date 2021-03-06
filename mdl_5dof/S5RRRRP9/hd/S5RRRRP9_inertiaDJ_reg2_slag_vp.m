% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP9_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:04
% EndTime: 2019-12-31 22:06:12
% DurationCPUTime: 2.53s
% Computational Cost: add. (2689->267), mult. (6690->458), div. (0->0), fcn. (5508->6), ass. (0->134)
t92 = cos(qJ(2));
t151 = t92 * qJD(2);
t89 = sin(qJ(3));
t136 = t89 * t151;
t91 = cos(qJ(3));
t154 = qJD(3) * t91;
t90 = sin(qJ(2));
t173 = t90 * t154 + t136;
t172 = -0.4e1 * t90;
t164 = cos(qJ(4));
t167 = t89 * pkin(6);
t132 = -pkin(3) - t167;
t160 = t90 * t91;
t165 = t92 * pkin(2);
t117 = -t90 * pkin(7) - t165;
t110 = -pkin(1) + t117;
t55 = t91 * t110;
t35 = -pkin(8) * t160 + t132 * t92 + t55;
t161 = t89 * t90;
t159 = t91 * t92;
t73 = pkin(6) * t159;
t45 = t89 * t110 + t73;
t39 = -pkin(8) * t161 + t45;
t88 = sin(qJ(4));
t171 = t164 * t39 + t88 * t35;
t125 = t164 * qJD(4);
t170 = t164 * qJD(3) + t125;
t84 = t89 ^ 2;
t86 = t91 ^ 2;
t157 = t84 - t86;
t124 = qJD(3) * t157;
t85 = t90 ^ 2;
t123 = (-t92 ^ 2 + t85) * qJD(2);
t169 = qJD(3) + qJD(4);
t153 = qJD(3) * t92;
t139 = t89 * t153;
t82 = t90 * qJD(2);
t101 = t91 * t82 + t139;
t116 = pkin(2) * t90 - pkin(7) * t92;
t108 = t116 * t89;
t29 = t101 * pkin(6) - qJD(2) * t108 - qJD(3) * t55;
t21 = -pkin(8) * t173 - t29;
t166 = t91 * pkin(2);
t168 = -pkin(8) - pkin(7);
t94 = (-t73 + (-t168 * t90 + pkin(1) + t165) * t89) * qJD(3) + (t168 * t159 + (-t132 + t166) * t90) * qJD(2);
t4 = -qJD(4) * t171 + t164 * t94 - t88 * t21;
t93 = 2 * qJD(5);
t163 = t88 * t89;
t162 = t88 * t91;
t59 = pkin(3) * t161 + t90 * pkin(6);
t155 = qJD(3) * t89;
t152 = qJD(4) * t88;
t150 = t92 * qJD(5);
t127 = qJD(2) * t164;
t118 = t92 * t127;
t140 = t90 * t155;
t144 = t88 * t161;
t23 = t89 * t118 - t88 * t140 - qJD(4) * t144 + (t88 * t151 + t170 * t90) * t91;
t57 = t164 * t89 + t162;
t46 = t57 * t90;
t149 = 0.2e1 * t46 * t23;
t38 = t169 * t57;
t131 = t164 * t91;
t56 = -t131 + t163;
t148 = 0.2e1 * t56 * t38;
t147 = t92 * t167;
t146 = -0.2e1 * pkin(1) * qJD(2);
t145 = -0.2e1 * pkin(2) * qJD(3);
t80 = pkin(6) * t151;
t43 = t173 * pkin(3) + t80;
t143 = pkin(4) * t82;
t142 = pkin(3) * t155;
t141 = pkin(3) * t152;
t137 = t91 * t153;
t135 = t89 * t154;
t134 = t90 * t151;
t133 = t91 * t151;
t79 = -t91 * pkin(3) - pkin(2);
t119 = t168 * t164;
t109 = qJD(3) * t119;
t111 = t89 * t119;
t129 = t168 * qJD(3);
t63 = t168 * t91;
t24 = -qJD(4) * t111 - t89 * t109 - t129 * t162 - t63 * t152;
t25 = -t63 * t125 - t91 * t109 + (qJD(4) * t168 + t129) * t163;
t40 = -t88 * t63 - t111;
t41 = t168 * t163 - t164 * t63;
t130 = -t41 * t24 + t40 * t25;
t122 = 0.2e1 * t134;
t121 = t89 * t133;
t120 = t85 * t135;
t22 = -t91 * t118 + t88 * t136 + t38 * t90;
t47 = t90 * t131 - t144;
t115 = t22 * t46 - t47 * t23;
t114 = t23 * t56 + t46 * t38;
t37 = t169 * t163 - t170 * t91;
t113 = t37 * t56 - t57 * t38;
t44 = t55 - t147;
t112 = -t44 * t91 - t45 * t89;
t107 = -t92 * t23 + t46 * t82;
t106 = t24 * t92 + t41 * t82;
t105 = t25 * t92 - t40 * t82;
t104 = -t92 * t38 + t56 * t82;
t18 = t164 * t35 - t88 * t39;
t3 = -t35 * t125 + t39 * t152 - t164 * t21 - t88 * t94;
t102 = -t40 * t22 - t41 * t23 + t24 * t46 + t25 * t47;
t76 = qJ(5) * t82;
t99 = -t3 + t76;
t98 = 0.2e1 * t24 * t56 + 0.2e1 * t25 * t57 - 0.2e1 * t40 * t37 - 0.2e1 * t41 * t38;
t97 = t22 * t56 - t57 * t23 + t37 * t46 - t47 * t38;
t30 = -t45 * qJD(3) + (pkin(6) * t161 + t91 * t116) * qJD(2);
t96 = t112 * qJD(3) - t29 * t91 - t30 * t89;
t95 = t92 * t141 + t4;
t81 = pkin(3) * t125;
t78 = -t164 * pkin(3) - pkin(4);
t75 = t88 * pkin(3) + qJ(5);
t74 = -0.2e1 * t141;
t69 = t81 + qJD(5);
t68 = -0.2e1 * t134;
t42 = t90 * t124 - t121;
t34 = t56 * pkin(4) - t57 * qJ(5) + t79;
t31 = -0.2e1 * t57 * t37;
t28 = t46 * pkin(4) - t47 * qJ(5) + t59;
t27 = t37 * t92 + t57 * t82;
t17 = t92 * pkin(4) - t18;
t16 = -t92 * qJ(5) + t171;
t11 = -0.2e1 * t47 * t22;
t10 = t38 * pkin(4) + t37 * qJ(5) - t57 * qJD(5) + t142;
t7 = 0.2e1 * t22 * t92 + 0.2e1 * t47 * t82;
t6 = -t22 * t57 - t47 * t37;
t5 = t23 * pkin(4) + t22 * qJ(5) - t47 * qJD(5) + t43;
t2 = -t143 - t4;
t1 = t99 - t150;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, -0.2e1 * t123, 0, t68, 0, 0, t90 * t146, t92 * t146, 0, 0, 0.2e1 * t86 * t134 - 0.2e1 * t120, t121 * t172 + 0.2e1 * t85 * t124, 0.2e1 * t91 * t123 + 0.2e1 * t90 * t139, 0.2e1 * t84 * t134 + 0.2e1 * t120, -0.2e1 * t89 * t123 + 0.2e1 * t90 * t137, t68, 0.2e1 * t44 * t82 - 0.2e1 * t30 * t92 + 0.2e1 * (t89 * t122 + t85 * t154) * pkin(6), -0.2e1 * t45 * t82 - 0.2e1 * t29 * t92 + 0.2e1 * (t91 * t122 - t85 * t155) * pkin(6), 0.2e1 * t112 * t151 + 0.2e1 * (t29 * t89 - t30 * t91 + (t44 * t89 - t45 * t91) * qJD(3)) * t90, 0.2e1 * pkin(6) ^ 2 * t134 - 0.2e1 * t45 * t29 + 0.2e1 * t44 * t30, t11, 0.2e1 * t115, t7, t149, -0.2e1 * t107, t68, 0.2e1 * t18 * t82 + 0.2e1 * t59 * t23 - 0.2e1 * t4 * t92 + 0.2e1 * t43 * t46, -0.2e1 * t171 * t82 - 0.2e1 * t59 * t22 - 0.2e1 * t3 * t92 + 0.2e1 * t43 * t47, -0.2e1 * t171 * t23 + 0.2e1 * t18 * t22 + 0.2e1 * t3 * t46 - 0.2e1 * t4 * t47, -0.2e1 * t171 * t3 + 0.2e1 * t18 * t4 + 0.2e1 * t59 * t43, t11, t7, -0.2e1 * t115, t68, 0.2e1 * t107, t149, -0.2e1 * t17 * t82 + 0.2e1 * t2 * t92 + 0.2e1 * t28 * t23 + 0.2e1 * t5 * t46, -0.2e1 * t1 * t46 - 0.2e1 * t16 * t23 - 0.2e1 * t17 * t22 + 0.2e1 * t2 * t47, -0.2e1 * t1 * t92 + 0.2e1 * t16 * t82 + 0.2e1 * t28 * t22 - 0.2e1 * t5 * t47, 0.2e1 * t16 * t1 + 0.2e1 * t17 * t2 + 0.2e1 * t28 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, 0, -t82, 0, -t80, pkin(6) * t82, 0, 0, -t42, t135 * t172 - t157 * t151, t89 * t82 - t137, t42, t101, 0, (pkin(7) * t159 + (-t166 + t167) * t90) * qJD(3) + (t117 * t89 - t73) * qJD(2), (pkin(6) * t160 + t108) * qJD(3) + (t117 * t91 + t147) * qJD(2), t96, -pkin(2) * t80 + t96 * pkin(7), t6, t97, t27, t114, -t104, 0, t46 * t142 + t79 * t23 + t59 * t38 + t43 * t56 + t105, t47 * t142 - t79 * t22 - t59 * t37 + t43 * t57 - t106, -t171 * t38 + t18 * t37 + t3 * t56 - t4 * t57 + t102, t59 * t142 - t171 * t24 - t18 * t25 - t3 * t41 - t4 * t40 + t43 * t79, t6, t27, -t97, 0, t104, t114, t10 * t46 + t34 * t23 + t28 * t38 + t5 * t56 + t105, -t1 * t56 - t16 * t38 - t17 * t37 + t2 * t57 + t102, -t10 * t47 + t34 * t22 + t28 * t37 - t5 * t57 + t106, t1 * t41 + t28 * t10 - t16 * t24 + t17 * t25 + t2 * t40 + t5 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t135, -0.2e1 * t124, 0, -0.2e1 * t135, 0, 0, t89 * t145, t91 * t145, 0, 0, t31, 0.2e1 * t113, 0, t148, 0, 0, 0.2e1 * t56 * t142 + 0.2e1 * t79 * t38, 0.2e1 * t57 * t142 - 0.2e1 * t79 * t37, t98, 0.2e1 * t79 * t142 + 0.2e1 * t130, t31, 0, -0.2e1 * t113, 0, 0, t148, 0.2e1 * t10 * t56 + 0.2e1 * t34 * t38, t98, -0.2e1 * t10 * t57 + 0.2e1 * t34 * t37, 0.2e1 * t34 * t10 + 0.2e1 * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133 - t140, 0, -t173, t82, t30, t29, 0, 0, 0, 0, -t22, 0, -t23, t82, pkin(3) * t90 * t127 + t95, (t92 * t125 - t88 * t82) * pkin(3) + t3, (t164 * t22 - t23 * t88 + (-t164 * t46 + t47 * t88) * qJD(4)) * pkin(3), (t164 * t4 - t3 * t88 + (t164 * t171 - t18 * t88) * qJD(4)) * pkin(3), 0, -t22, 0, t82, t23, 0, (pkin(4) - t78) * t82 + t95, t47 * t141 - t78 * t22 - t75 * t23 - t69 * t46, t75 * t82 + (-qJD(5) - t69) * t92 + t99, t1 * t75 + t141 * t17 + t16 * t69 + t2 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, 0, -t155, 0, -pkin(7) * t154, pkin(7) * t155, 0, 0, 0, 0, -t37, 0, -t38, 0, -t25, t24, (t164 * t37 - t38 * t88 + (-t164 * t56 + t57 * t88) * qJD(4)) * pkin(3), (-t164 * t25 - t24 * t88 + (t164 * t41 + t40 * t88) * qJD(4)) * pkin(3), 0, -t37, 0, 0, t38, 0, -t25, t57 * t141 - t78 * t37 - t75 * t38 - t69 * t56, -t24, t141 * t40 - t24 * t75 + t25 * t78 + t41 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -0.2e1 * t81, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0.2e1 * t69, 0.2e1 * t141 * t78 + 0.2e1 * t75 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, -t23, t82, t4, t3, 0, 0, 0, -t22, 0, t82, t23, 0, t4 + 0.2e1 * t143, pkin(4) * t22 - t23 * qJ(5) - t46 * qJD(5), -t3 + 0.2e1 * t76 - 0.2e1 * t150, -t2 * pkin(4) + t1 * qJ(5) + qJD(5) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, -t38, 0, -t25, t24, 0, 0, 0, -t37, 0, 0, t38, 0, -t25, pkin(4) * t37 - t38 * qJ(5) - t56 * qJD(5), -t24, -t25 * pkin(4) - t24 * qJ(5) + t41 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, -t81, 0, 0, 0, 0, 0, 0, 0, 0, -t141, 0, t93 + t81, -pkin(4) * t141 + t69 * qJ(5) + t75 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, qJ(5) * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t22, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
