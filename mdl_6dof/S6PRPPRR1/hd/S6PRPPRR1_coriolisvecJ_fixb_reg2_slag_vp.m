% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPPRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:15
% EndTime: 2019-03-08 19:16:20
% DurationCPUTime: 2.23s
% Computational Cost: add. (4432->285), mult. (11196->395), div. (0->0), fcn. (9151->12), ass. (0->166)
t114 = sin(pkin(11));
t121 = sin(qJ(2));
t115 = sin(pkin(6));
t169 = qJD(1) * t115;
t158 = t121 * t169;
t102 = t114 * t158;
t117 = cos(pkin(11));
t123 = cos(qJ(2));
t157 = t123 * t169;
t75 = t117 * t157 - t102;
t208 = t75 - qJD(4);
t113 = sin(pkin(12));
t120 = sin(qJ(5));
t116 = cos(pkin(12));
t201 = cos(qJ(5));
t159 = t201 * t116;
t127 = -t120 * t113 + t159;
t107 = t114 * pkin(2) + qJ(4);
t187 = pkin(8) + t107;
t90 = t187 * t113;
t91 = t187 * t116;
t129 = -t120 * t91 - t201 * t90;
t185 = -t129 * qJD(5) + t208 * t127;
t77 = (t114 * t123 + t117 * t121) * t115;
t72 = qJD(1) * t77;
t88 = t127 * qJD(5);
t94 = t201 * t113 + t120 * t116;
t89 = t94 * qJD(5);
t213 = -t89 * pkin(5) + t88 * pkin(9) + t72;
t167 = qJD(2) * t115;
t155 = qJD(1) * t167;
t145 = t123 * t155;
t97 = t117 * t145;
t62 = t97 + (qJD(4) - t102) * qJD(2);
t212 = t94 * t62;
t119 = sin(qJ(6));
t122 = cos(qJ(6));
t118 = cos(pkin(6));
t106 = t118 * qJD(1) + qJD(3);
t100 = t116 * t106;
t98 = qJD(2) * pkin(2) + t157;
t64 = t114 * t98 + t117 * t158;
t60 = qJD(2) * qJ(4) + t64;
t37 = t100 + (-pkin(8) * qJD(2) - t60) * t113;
t166 = qJD(2) * t116;
t42 = t113 * t106 + t116 * t60;
t38 = pkin(8) * t166 + t42;
t18 = t120 * t37 + t201 * t38;
t16 = qJD(5) * pkin(9) + t18;
t63 = t117 * t98 - t102;
t137 = qJD(4) - t63;
t160 = -t116 * pkin(4) - pkin(3);
t50 = t160 * qJD(2) + t137;
t104 = qJD(2) * t159;
t168 = qJD(2) * t113;
t156 = t120 * t168;
t84 = -t104 + t156;
t86 = t94 * qJD(2);
t25 = t84 * pkin(5) - t86 * pkin(9) + t50;
t134 = t119 * t16 - t122 * t25;
t146 = t121 * t155;
t65 = t114 * t145 + t117 * t146;
t103 = qJD(5) * t104;
t78 = qJD(5) * t156 - t103;
t79 = qJD(2) * t89;
t30 = t79 * pkin(5) + t78 * pkin(9) + t65;
t131 = t120 * t38 - t201 * t37;
t210 = t127 * t62;
t9 = -t131 * qJD(5) + t210;
t1 = -t134 * qJD(6) + t119 * t30 + t122 * t9;
t82 = qJD(6) + t84;
t144 = t134 * t82 + t1;
t6 = t119 * t25 + t122 * t16;
t2 = -qJD(6) * t6 - t119 * t9 + t122 * t30;
t211 = t6 * t82 + t2;
t151 = t119 * t82;
t69 = t119 * qJD(5) + t122 * t86;
t209 = t69 * t151;
t153 = t72 * qJD(2) - t65;
t176 = qJD(6) * t69;
t40 = -t119 * t78 + t176;
t56 = t94 * t79;
t138 = t82 * t88 + t56;
t165 = qJD(6) * t119;
t162 = t94 * t165;
t207 = -t122 * t138 + t82 * t162;
t206 = t86 ^ 2;
t196 = t117 * pkin(2);
t101 = t160 - t196;
t46 = -pkin(5) * t127 - t94 * pkin(9) + t101;
t48 = -t120 * t90 + t201 * t91;
t24 = t119 * t46 + t122 * t48;
t203 = t24 * qJD(6) - t185 * t119 + t213 * t122;
t23 = -t119 * t48 + t122 * t46;
t202 = -t23 * qJD(6) + t213 * t119 + t185 * t122;
t10 = t18 * qJD(5) + t212;
t58 = -t77 * t113 + t118 * t116;
t59 = t118 * t113 + t77 * t116;
t130 = -t120 * t59 + t201 * t58;
t200 = t10 * t130;
t199 = t10 * t129;
t198 = t10 * t127;
t197 = t10 * t94;
t133 = t114 * t121 - t117 * t123;
t76 = t133 * t115;
t43 = t65 * t76;
t163 = t122 * qJD(5);
t67 = t119 * t86 - t163;
t195 = t67 * t84;
t194 = t67 * t88;
t193 = t69 * t67;
t192 = t69 * t86;
t191 = t69 * t88;
t190 = t79 * t127;
t189 = t86 * t67;
t188 = t86 * t84;
t186 = -t48 * qJD(5) + t208 * t94;
t164 = qJD(6) * t122;
t184 = -t119 * t40 - t67 * t164;
t178 = t40 * t122;
t183 = -t122 * t194 - t94 * t178;
t39 = -qJD(6) * t163 + t122 * t78 + t86 * t165;
t182 = t127 * t39 + t69 * t89;
t181 = -t88 * t84 - t56;
t179 = t119 * t79;
t177 = qJD(6) * t67;
t124 = qJD(2) ^ 2;
t175 = t115 * t124;
t74 = t133 * t167;
t172 = t74 * qJD(2);
t171 = t88 * qJD(5);
t170 = t113 ^ 2 + t116 ^ 2;
t161 = t94 * t164;
t154 = t170 * t62;
t152 = -t39 + t177;
t150 = t122 * t82;
t148 = t69 * t161;
t147 = pkin(9) * qJD(6) * t82 + t10;
t15 = -qJD(5) * pkin(5) + t131;
t143 = t15 * t88 + t197;
t142 = t119 * t6 - t122 * t134;
t141 = -t119 * t134 - t122 * t6;
t140 = t127 * t40 - t89 * t67;
t139 = t127 * t78 + t86 * t89;
t41 = -t113 * t60 + t100;
t135 = t41 * t113 - t42 * t116;
t28 = t120 * t58 + t201 * t59;
t20 = t76 * t119 + t122 * t28;
t19 = -t119 * t28 + t76 * t122;
t132 = t122 * t79 - t84 * t151 - t82 * t165;
t128 = -pkin(9) * t79 + t82 * t15;
t126 = -t138 * t119 - t82 * t161;
t125 = -t142 * qJD(6) + t1 * t122 - t2 * t119;
t83 = t84 ^ 2;
t81 = t89 * qJD(5);
t73 = qJD(2) * t77;
t66 = -t114 * t146 + t97;
t57 = -qJD(2) * pkin(3) + t137;
t51 = t86 * pkin(5) + t84 * pkin(9);
t14 = t28 * qJD(5) - t94 * t74;
t13 = t130 * qJD(5) - t127 * t74;
t12 = t119 * t51 - t122 * t131;
t11 = t119 * t131 + t122 * t51;
t4 = -t20 * qJD(6) - t119 * t13 + t73 * t122;
t3 = t19 * qJD(6) + t73 * t119 + t122 * t13;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121 * t175, -t123 * t175, 0, 0, 0, 0, 0, 0, 0, 0, -t73 * qJD(2), t172, 0, -t63 * t73 - t64 * t74 + t66 * t77 + t43, 0, 0, 0, 0, 0, 0, -t73 * t166, t73 * t168, -t170 * t172, t57 * t73 + t43 + (-t42 * t74 + t59 * t62) * t116 + (t41 * t74 - t58 * t62) * t113, 0, 0, 0, 0, 0, 0, -t14 * qJD(5) + t73 * t84 + t76 * t79, -t13 * qJD(5) + t73 * t86 - t76 * t78, -t13 * t84 + t130 * t78 + t14 * t86 - t28 * t79, t18 * t13 + t131 * t14 + t9 * t28 + t50 * t73 - t200 + t43, 0, 0, 0, 0, 0, 0, -t130 * t40 + t14 * t67 + t19 * t79 + t4 * t82, t130 * t39 + t14 * t69 - t20 * t79 - t3 * t82, t19 * t39 - t20 * t40 - t3 * t67 - t4 * t69, t1 * t20 - t134 * t4 + t15 * t14 + t2 * t19 + t6 * t3 - t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, -t97 + (t75 + t102) * qJD(2), 0, t63 * t72 - t64 * t75 + (t114 * t66 - t117 * t65) * pkin(2), 0, 0, 0, 0, 0, 0, t153 * t116, -t153 * t113, -t208 * qJD(2) * t170 + t154, t65 * (-pkin(3) - t196) - t57 * t72 + t107 * t154 + t208 * t135, -t78 * t94 + t86 * t88, -t139 + t181, t171, t84 * t89 - t190, -t81, 0, t186 * qJD(5) + t101 * t79 - t127 * t65 + t50 * t89 - t72 * t84, t185 * qJD(5) - t101 * t78 + t50 * t88 + t65 * t94 - t72 * t86, t127 * t9 + t129 * t78 + t131 * t88 - t18 * t89 + t185 * t84 - t186 * t86 - t48 * t79 + t197, t65 * t101 - t131 * t186 - t185 * t18 + t9 * t48 - t50 * t72 - t199, -t69 * t162 + (-t39 * t94 + t191) * t122, -t148 + (-t191 + (t39 + t177) * t94) * t119 + t183, t182 - t207, t67 * t161 + (t40 * t94 + t194) * t119, t126 + t140, t82 * t89 - t190, t143 * t119 - t127 * t2 - t129 * t40 - t134 * t89 + t15 * t161 - t186 * t67 - t203 * t82 + t23 * t79, t1 * t127 + t143 * t122 + t129 * t39 - t15 * t162 - t186 * t69 + t202 * t82 - t24 * t79 - t6 * t89, t23 * t39 - t24 * t40 - t142 * t88 + t203 * t69 + t202 * t67 + (qJD(6) * t141 - t1 * t119 - t2 * t122) * t94, t1 * t24 + t134 * t203 - t186 * t15 + t2 * t23 - t202 * t6 - t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t171, t139 + t181, t131 * t89 + t18 * t88 + t9 * t94 - t198, 0, 0, 0, 0, 0, 0, t126 - t140, t182 + t207, t148 + (t152 * t94 + t191) * t119 + t183, t125 * t94 - t141 * t88 + t15 * t89 - t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170 * t124, t135 * qJD(2) + t65, 0, 0, 0, 0, 0, 0, 0.2e1 * t86 * qJD(5), t103 + (-t84 - t156) * qJD(5), -t83 - t206, -t131 * t86 + t18 * t84 + t65, 0, 0, 0, 0, 0, 0, t132 - t189, -t82 ^ 2 * t122 - t179 - t192 (t39 - t195) * t122 + t209 + t184, t144 * t119 + t211 * t122 - t15 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t188, -t83 + t206, t103 + (t84 - t156) * qJD(5), -t188, 0, 0, -t50 * t86 - t212, t50 * t84 - t210, 0, 0, -t39 * t119 + t150 * t69 (-t39 - t195) * t122 - t209 + t184, t150 * t82 + t179 - t192, t151 * t67 - t178, t132 + t189, -t82 * t86, -pkin(5) * t40 - t11 * t82 + t119 * t128 - t122 * t147 + t134 * t86 - t18 * t67, pkin(5) * t39 + t119 * t147 + t12 * t82 + t122 * t128 - t18 * t69 + t6 * t86, t11 * t69 + t12 * t67 + ((-t40 + t176) * pkin(9) + t144) * t122 + (t152 * pkin(9) - t211) * t119, -t10 * pkin(5) + pkin(9) * t125 + t11 * t134 - t6 * t12 - t15 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t67 ^ 2 + t69 ^ 2, t67 * t82 - t39, -t193, t69 * t82 - t40, t79, -t15 * t69 + t211, t15 * t67 - t144, 0, 0;];
tauc_reg  = t5;
