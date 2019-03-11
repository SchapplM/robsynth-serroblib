% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPPRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:01
% EndTime: 2019-03-08 19:20:07
% DurationCPUTime: 2.08s
% Computational Cost: add. (2477->257), mult. (5797->370), div. (0->0), fcn. (4298->10), ass. (0->154)
t73 = sin(pkin(11));
t74 = sin(pkin(6));
t75 = cos(pkin(11));
t79 = sin(qJ(2));
t82 = cos(qJ(2));
t49 = (-t73 * t79 + t75 * t82) * t74;
t147 = qJD(1) * t74;
t126 = t82 * t147;
t61 = qJD(2) * pkin(2) + t126;
t127 = t79 * t147;
t62 = t73 * t127;
t37 = t75 * t61 - t62;
t102 = qJD(4) - t37;
t26 = (-pkin(3) - pkin(8)) * qJD(2) + t102;
t76 = cos(pkin(6));
t65 = qJD(1) * t76 + qJD(3);
t78 = sin(qJ(5));
t81 = cos(qJ(5));
t22 = t26 * t78 + t65 * t81;
t20 = qJD(5) * pkin(9) + t22;
t38 = t75 * t127 + t73 * t61;
t92 = pkin(5) * t78 - pkin(9) * t81 + qJ(4);
t25 = t92 * qJD(2) + t38;
t77 = sin(qJ(6));
t80 = cos(qJ(6));
t100 = t20 * t77 - t25 * t80;
t47 = qJD(2) * t49;
t41 = t47 * qJD(1);
t107 = pkin(5) * t81 + pkin(9) * t78;
t54 = t107 * qJD(5) + qJD(4);
t24 = t54 * qJD(2) + t41;
t21 = t26 * t81 - t65 * t78;
t95 = t73 * t82 + t75 * t79;
t50 = t95 * t74;
t46 = qJD(2) * t50;
t40 = qJD(1) * t46;
t7 = t21 * qJD(5) + t78 * t40;
t1 = -t100 * qJD(6) + t77 * t24 + t80 * t7;
t144 = qJD(2) * t78;
t67 = qJD(6) + t144;
t185 = t100 * t67 + t1;
t48 = t75 * t126 - t62;
t184 = -t48 + t54;
t141 = qJD(5) * t78;
t120 = t77 * t141;
t135 = qJD(6) * t81;
t183 = t80 * t135 - t120;
t6 = t20 * t80 + t25 * t77;
t2 = -qJD(6) * t6 + t80 * t24 - t77 * t7;
t182 = -t6 * t67 - t2;
t31 = qJD(2) * qJ(4) + t38;
t45 = qJD(1) * t50;
t181 = -t31 + t45;
t140 = qJD(5) * t81;
t131 = t80 * qJD(5);
t118 = t77 * t135;
t89 = t78 * t131 + t118;
t42 = t89 * qJD(2) - qJD(6) * t131;
t143 = qJD(2) * t81;
t124 = t80 * t143;
t142 = qJD(5) * t77;
t57 = t124 + t142;
t150 = t57 * t140 - t42 * t78;
t138 = qJD(6) * t57;
t43 = -qJD(2) * t120 + t138;
t105 = -t100 * t77 - t6 * t80;
t35 = t81 * t40;
t8 = t22 * qJD(5) - t35;
t180 = qJD(5) * t105 + t8;
t72 = t81 ^ 2;
t146 = qJD(2) * t72;
t179 = -(-t67 * t78 + t146) * t131 + t67 * t118;
t172 = t8 * t81;
t86 = -(t21 * t78 - t22 * t81) * qJD(5) + t7 * t78 - t172;
t178 = pkin(2) * t73;
t96 = -t49 * t81 - t76 * t78;
t177 = t96 * t8;
t174 = t77 * t8;
t173 = t8 * t80;
t19 = -qJD(5) * pkin(5) - t21;
t171 = t19 * t77;
t170 = t19 * t80;
t169 = t40 * t49;
t55 = t77 * t143 - t131;
t168 = t55 * t67;
t167 = t55 * t81;
t166 = t57 * t55;
t165 = t57 * t67;
t164 = t67 * t77;
t163 = t67 * t80;
t84 = qJD(2) ^ 2;
t161 = t74 * t84;
t159 = t77 * t78;
t158 = t77 * t81;
t157 = t78 * t43;
t156 = t78 * t80;
t155 = t80 * t81;
t83 = qJD(5) ^ 2;
t154 = t83 * t78;
t153 = t83 * t81;
t116 = -pkin(2) * t75 - pkin(3);
t66 = -pkin(8) + t116;
t121 = t66 * t140;
t52 = t92 + t178;
t27 = -t66 * t159 + t52 * t80;
t152 = t27 * qJD(6) + t80 * t121 - t45 * t156 + t184 * t77;
t28 = t66 * t156 + t52 * t77;
t151 = -t28 * qJD(6) - t77 * t121 + t45 * t159 + t184 * t80;
t71 = t78 ^ 2;
t149 = t71 - t72;
t148 = -t83 - t84;
t139 = qJD(6) * t55;
t137 = qJD(6) * t77;
t136 = qJD(6) * t80;
t134 = t31 * qJD(2);
t133 = t46 * qJD(2);
t132 = t47 * qJD(2);
t130 = qJD(4) - t48;
t129 = qJD(2) * qJD(5);
t128 = t81 * t84 * t78;
t125 = t77 * t146;
t123 = t55 * t141;
t122 = t57 * t141;
t119 = t57 * t135;
t115 = t81 * t129;
t114 = t67 + t144;
t113 = -t42 + t139;
t109 = qJD(6) * t78 + qJD(2);
t108 = t78 * t115;
t106 = -t100 * t80 + t6 * t77;
t68 = qJ(4) + t178;
t103 = qJD(2) * t68 - t181;
t99 = -t21 * t81 - t22 * t78;
t34 = qJD(4) * qJD(2) + t41;
t97 = t31 * t47 + t34 * t50;
t33 = -t49 * t78 + t76 * t81;
t16 = t33 * t80 + t50 * t77;
t15 = -t33 * t77 + t50 * t80;
t93 = t183 * t67;
t91 = -pkin(9) * t140 + t19 * t78;
t90 = t130 * t31 + t34 * t68;
t88 = t130 * qJD(2) - t66 * t83 + t34;
t87 = -t106 * qJD(6) + t1 * t80 - t2 * t77;
t85 = qJD(5) * t19 + t87;
t60 = t107 * qJD(2);
t30 = t43 * t155;
t29 = -qJD(2) * pkin(3) + t102;
t14 = t33 * qJD(5) - t46 * t81;
t13 = t96 * qJD(5) + t46 * t78;
t10 = t21 * t80 + t60 * t77;
t9 = -t21 * t77 + t60 * t80;
t4 = t15 * qJD(6) + t13 * t80 + t47 * t77;
t3 = -t16 * qJD(6) - t13 * t77 + t47 * t80;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79 * t161, -t82 * t161, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t132, 0, -t37 * t46 + t38 * t47 + t41 * t50 - t169, 0, 0, 0, 0, 0, 0, 0, t133, t132, t29 * t46 - t169 + t97, 0, 0, 0, 0, 0, 0, -t14 * qJD(5) + (t50 * t140 + t47 * t78) * qJD(2), -t13 * qJD(5) + (-t50 * t141 + t47 * t81) * qJD(2) (-t13 * t78 + t14 * t81 + (-t33 * t81 + t78 * t96) * qJD(5)) * qJD(2), t13 * t22 - t14 * t21 + t33 * t7 - t177 + t97, 0, 0, 0, 0, 0, 0, t115 * t15 + t14 * t55 + t3 * t67 - t43 * t96, -t115 * t16 + t14 * t57 - t4 * t67 + t42 * t96, t15 * t42 - t16 * t43 - t3 * t57 - t4 * t55, t1 * t16 - t100 * t3 + t14 * t19 + t15 * t2 + t4 * t6 - t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t95 * t147 + t45) * qJD(2), qJD(2) * t48 - t41, 0, t37 * t45 - t38 * t48 + (-t40 * t75 + t41 * t73) * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0 (0.2e1 * qJD(4) - t48) * qJD(2) + t41, t40 * t116 - t29 * t45 + t90, -0.2e1 * t108, 0.2e1 * t149 * t129, -t154, 0.2e1 * t108, -t153, 0, t103 * t140 + t88 * t78, -t103 * t141 + t88 * t81 (t71 + t72) * t45 * qJD(2) - t86, t99 * t45 + t86 * t66 + t90, -t42 * t155 - t89 * t57, -t30 + (-t119 + t123) * t80 + (t122 + (t42 + t139) * t81) * t77, t150 - t179, t43 * t158 + t183 * t55, -t157 + (-t125 - t167) * qJD(5) - t93, t114 * t140, t151 * t67 + (t2 + (t55 * t66 - t171) * qJD(5)) * t78 + (t19 * t136 - t43 * t66 + t45 * t55 + t174 + (qJD(2) * t27 - t100) * qJD(5)) * t81, -t152 * t67 + (-t1 + (t57 * t66 - t170) * qJD(5)) * t78 + (-t19 * t137 + t42 * t66 + t45 * t57 + t173 + (-qJD(2) * t28 - t6) * qJD(5)) * t81, t27 * t42 - t28 * t43 - t151 * t57 - t152 * t55 + t106 * t141 + (qJD(6) * t105 - t1 * t77 - t2 * t80) * t81, -t66 * t172 + t1 * t28 + t2 * t27 + t152 * t6 - t151 * t100 + (t141 * t66 + t45 * t81) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t154, 0, t99 * qJD(5) + t7 * t81 + t8 * t78, 0, 0, 0, 0, 0, 0, t157 + (-t125 + t167) * qJD(5) - t93, t150 + t179, -t30 + (t119 + t123) * t80 + (t113 * t81 - t122) * t77, t180 * t78 + t85 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t181 * qJD(2), 0, 0, 0, 0, 0, 0, t148 * t78, t148 * t81, 0, t86 - t134, 0, 0, 0, 0, 0, 0, -t43 * t81 - t109 * t163 + (-t114 * t158 + t55 * t78) * qJD(5), t42 * t81 + t109 * t164 + (-t67 * t155 + (t57 - t124) * t78) * qJD(5) (t109 * t57 - t140 * t55 - t157) * t80 + (t109 * t55 + t150) * t77, -t106 * qJD(2) - t180 * t81 + t85 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, -t149 * t84, 0, -t128, 0, 0, -t81 * t134 + t35 (-t40 + t134) * t78, 0, 0, t57 * t163 - t42 * t77 (-t42 - t168) * t80 + (-t43 - t165) * t77, t67 * t136 + (t67 * t156 + (-t57 + t142) * t81) * qJD(2), t164 * t55 - t43 * t80, -t67 * t137 + (-t67 * t159 + (t55 + t131) * t81) * qJD(2), -t67 * t143, -pkin(5) * t43 - t22 * t55 - t67 * t9 - t173 + (-pkin(9) * t163 + t171) * qJD(6) + (t100 * t81 + t77 * t91) * qJD(2), pkin(5) * t42 + t10 * t67 - t22 * t57 + t174 + (pkin(9) * t164 + t170) * qJD(6) + (t6 * t81 + t80 * t91) * qJD(2), t10 * t55 + t57 * t9 + ((-t43 + t138) * pkin(9) + t185) * t80 + (pkin(9) * t113 + t182) * t77, -pkin(5) * t8 + pkin(9) * t87 - t10 * t6 + t100 * t9 - t19 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, -t55 ^ 2 + t57 ^ 2, t168 - t42, -t166, t165 - t43, t115, -t19 * t57 - t182, t19 * t55 - t185, 0, 0;];
tauc_reg  = t5;
