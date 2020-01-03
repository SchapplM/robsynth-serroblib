% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP9
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
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:28
% EndTime: 2019-12-31 20:07:35
% DurationCPUTime: 2.03s
% Computational Cost: add. (2774->296), mult. (7014->408), div. (0->0), fcn. (4769->6), ass. (0->157)
t138 = cos(qJ(2));
t186 = qJD(1) * t138;
t123 = -qJD(4) + t186;
t135 = sin(qJ(4));
t137 = cos(qJ(4));
t133 = sin(pkin(8));
t136 = sin(qJ(2));
t187 = qJD(1) * t136;
t171 = t133 * t187;
t134 = cos(pkin(8));
t184 = qJD(2) * t134;
t96 = -t171 + t184;
t175 = t134 * t187;
t185 = qJD(2) * t133;
t97 = t175 + t185;
t51 = t135 * t97 - t137 * t96;
t201 = t123 * t51;
t179 = qJD(1) * qJD(2);
t170 = t138 * t179;
t160 = t133 * t170;
t182 = qJD(2) * t138;
t192 = t137 * t134;
t162 = t182 * t192;
t180 = qJD(4) * t137;
t21 = (qJD(4) * t97 + t160) * t135 - qJD(1) * t162 - t96 * t180;
t221 = t21 - t201;
t220 = t51 ^ 2;
t181 = qJD(4) * t135;
t215 = -t133 * t181 + t134 * t180;
t99 = t133 * t135 - t192;
t204 = t99 * t186 + t215;
t100 = t133 * t137 + t134 * t135;
t146 = t100 * t138;
t87 = t100 * qJD(4);
t203 = -qJD(1) * t146 + t87;
t153 = t135 * t96 + t137 * t97;
t212 = t153 ^ 2;
t219 = -0.2e1 * t179;
t193 = t134 * t138;
t151 = pkin(3) * t136 - pkin(7) * t193;
t158 = pkin(2) * t136 - qJ(3) * t138;
t102 = t158 * qJD(1);
t67 = pkin(6) * t171 + t134 * t102;
t44 = t151 * qJD(1) + t67;
t194 = t134 * t136;
t195 = t133 * t138;
t147 = -pkin(6) * t194 - pkin(7) * t195;
t88 = t133 * t102;
t56 = t147 * qJD(1) + t88;
t209 = pkin(7) + qJ(3);
t110 = t209 * t133;
t111 = t209 * t134;
t65 = -t110 * t135 + t111 * t137;
t218 = t100 * qJD(3) + t65 * qJD(4) - t135 * t56 + t137 * t44;
t152 = -t110 * t137 - t111 * t135;
t217 = -t99 * qJD(3) + t152 * qJD(4) - t135 * t44 - t137 * t56;
t216 = t123 * t153;
t214 = t153 * qJD(4);
t107 = -pkin(2) * t138 - qJ(3) * t136 - pkin(1);
t95 = t134 * t107;
t57 = -pkin(7) * t194 + t95 + (-pkin(6) * t133 - pkin(3)) * t138;
t196 = t133 * t136;
t121 = pkin(6) * t193;
t72 = t133 * t107 + t121;
t63 = -pkin(7) * t196 + t72;
t154 = t135 * t57 + t137 * t63;
t145 = t151 * qJD(2);
t183 = qJD(2) * t136;
t176 = pkin(6) * t183;
t84 = t158 * qJD(2) - t136 * qJD(3);
t61 = t133 * t176 + t134 * t84;
t37 = t145 + t61;
t76 = t133 * t84;
t45 = t147 * qJD(2) + t76;
t213 = -t154 * qJD(4) - t135 * t45 + t137 * t37;
t127 = pkin(6) * t187;
t202 = qJD(2) * pkin(2);
t166 = qJD(3) - t202;
t106 = t127 + t166;
t66 = -pkin(3) * t96 + t106;
t12 = pkin(4) * t51 - qJ(5) * t153 + t66;
t211 = t12 * t153;
t210 = t153 * t51;
t208 = qJ(5) * t187 - t217;
t207 = pkin(4) * t187 + t218;
t128 = pkin(6) * t186;
t177 = pkin(3) * t186;
t91 = t133 * t177 + t128;
t206 = -t203 * pkin(4) + t204 * qJ(5) + qJD(5) * t100 + t91;
t105 = (qJD(3) - t127) * qJD(2);
t75 = t84 * qJD(1);
t43 = t134 * t105 + t133 * t75;
t113 = qJD(2) * qJ(3) + t128;
t90 = t107 * qJD(1);
t58 = -t113 * t133 + t134 * t90;
t28 = -pkin(7) * t97 - t177 + t58;
t59 = t134 * t113 + t133 * t90;
t33 = pkin(7) * t96 + t59;
t9 = -t135 * t33 + t137 * t28;
t200 = qJD(5) - t9;
t198 = qJD(2) * t152;
t197 = qJD(2) * t65;
t140 = qJD(1) ^ 2;
t191 = t138 * t140;
t139 = qJD(2) ^ 2;
t190 = t139 * t136;
t189 = t139 * t138;
t122 = pkin(6) * t170;
t83 = pkin(3) * t160 + t122;
t129 = pkin(6) * t182;
t174 = t133 * t182;
t92 = pkin(3) * t174 + t129;
t103 = pkin(3) * t196 + t136 * pkin(6);
t188 = t136 ^ 2 - t138 ^ 2;
t178 = pkin(6) * t195;
t126 = -pkin(3) * t134 - pkin(2);
t169 = t136 * t179;
t42 = -t105 * t133 + t134 * t75;
t25 = qJD(1) * t145 + t42;
t32 = -pkin(7) * t160 + t43;
t168 = t135 * t32 - t137 * t25 + t33 * t180 + t28 * t181;
t167 = pkin(1) * t219;
t165 = -t97 + t185;
t164 = -t96 + t184;
t163 = pkin(4) * t169;
t161 = qJ(5) * t169;
t159 = -t106 + t166;
t10 = t135 * t28 + t137 * t33;
t155 = -t135 * t63 + t137 * t57;
t150 = -t10 * t123 - t168;
t149 = -t135 * t25 - t137 * t32 - t28 * t180 + t33 * t181;
t148 = t135 * t37 + t137 * t45 + t57 * t180 - t63 * t181;
t143 = qJD(2) * t146;
t2 = -t163 + t168;
t22 = qJD(1) * t143 + t214;
t3 = pkin(4) * t22 + qJ(5) * t21 - qJD(5) * t153 + t83;
t141 = t22 - t216;
t81 = t99 * t136;
t80 = t100 * t136;
t71 = t95 - t178;
t68 = -pkin(6) * t175 + t88;
t62 = -t134 * t176 + t76;
t49 = pkin(4) * t99 - qJ(5) * t100 + t126;
t40 = t215 * t136 + t143;
t39 = t135 * t174 + t136 * t87 - t162;
t29 = pkin(4) * t80 + qJ(5) * t81 + t103;
t17 = pkin(4) * t153 + qJ(5) * t51;
t16 = pkin(4) * t138 - t155;
t15 = -qJ(5) * t138 + t154;
t11 = -t21 - t201;
t8 = -qJ(5) * t123 + t10;
t7 = pkin(4) * t123 + t200;
t6 = pkin(4) * t40 + qJ(5) * t39 + qJD(5) * t81 + t92;
t5 = -pkin(4) * t183 - t213;
t4 = qJ(5) * t183 - qJD(5) * t138 + t148;
t1 = -qJD(5) * t123 - t149 + t161;
t13 = [0, 0, 0, 0.2e1 * t138 * t169, t188 * t219, t189, -t190, 0, -pkin(6) * t189 + t136 * t167, pkin(6) * t190 + t138 * t167, (-qJD(1) * t61 - t42) * t138 + ((-pkin(6) * t96 + t106 * t133) * t138 + (t58 + (t71 + 0.2e1 * t178) * qJD(1)) * t136) * qJD(2), (qJD(1) * t62 + t43) * t138 + ((pkin(6) * t97 + t106 * t134) * t138 + (-t59 + (-t72 + 0.2e1 * t121) * qJD(1)) * t136) * qJD(2), -t61 * t97 + t62 * t96 + (-t133 * t43 - t134 * t42) * t136 + (-t133 * t59 - t134 * t58 + (-t133 * t72 - t134 * t71) * qJD(1)) * t182, t42 * t71 + t43 * t72 + t58 * t61 + t59 * t62 + (t106 + t127) * t129, -t153 * t39 + t21 * t81, -t153 * t40 + t21 * t80 + t22 * t81 + t39 * t51, t39 * t123 + t21 * t138 + (-qJD(1) * t81 + t153) * t183, t40 * t123 + t22 * t138 + (-qJD(1) * t80 - t51) * t183, (-t123 - t186) * t183, -t213 * t123 + t168 * t138 + t92 * t51 + t103 * t22 + t83 * t80 + t66 * t40 + (t155 * qJD(1) + t9) * t183, t148 * t123 - t149 * t138 + t92 * t153 - t103 * t21 - t83 * t81 - t66 * t39 + (-t154 * qJD(1) - t10) * t183, t12 * t40 + t123 * t5 + t138 * t2 + t22 * t29 + t3 * t80 + t51 * t6 + (-qJD(1) * t16 - t7) * t183, -t1 * t80 - t15 * t22 + t153 * t5 - t16 * t21 - t2 * t81 - t39 * t7 - t4 * t51 - t40 * t8, -t1 * t138 + t12 * t39 - t123 * t4 + t21 * t29 + t3 * t81 - t153 * t6 + (qJD(1) * t15 + t8) * t183, t1 * t15 + t12 * t6 + t16 * t2 + t29 * t3 + t4 * t8 + t5 * t7; 0, 0, 0, -t136 * t191, t188 * t140, 0, 0, 0, t140 * pkin(1) * t136, pkin(1) * t191, ((-qJ(3) * t185 - t58) * t136 + (-t164 * pkin(6) + t159 * t133 + t67) * t138) * qJD(1), ((-qJ(3) * t184 + t59) * t136 + (t165 * pkin(6) + t159 * t134 - t68) * t138) * qJD(1), t67 * t97 - t68 * t96 + (qJD(3) * t96 + t58 * t186 + t43) * t134 + (qJD(3) * t97 + t59 * t186 - t42) * t133, -t58 * t67 - t59 * t68 + (-t133 * t58 + t134 * t59) * qJD(3) + (-t133 * t42 + t134 * t43) * qJ(3) + (-t106 - t202) * t128, -t21 * t100 + t153 * t204, -t100 * t22 - t153 * t203 - t204 * t51 + t21 * t99, -t204 * t123 + (qJD(2) * t100 - t153) * t187, t203 * t123 + (-qJD(2) * t99 + t51) * t187, t123 * t187, t126 * t22 - t91 * t51 + t83 * t99 + t203 * t66 + t218 * t123 + (-t9 + t198) * t187, t83 * t100 - t126 * t21 - t91 * t153 + t204 * t66 + t217 * t123 + (t10 - t197) * t187, t22 * t49 + t3 * t99 - t206 * t51 + t207 * t123 + t203 * t12 + (t7 + t198) * t187, -t1 * t99 + t100 * t2 + t152 * t21 + t153 * t207 - t203 * t8 + t204 * t7 + t208 * t51 - t22 * t65, -t100 * t3 + t21 * t49 + t206 * t153 + t208 * t123 - t204 * t12 + (-t8 + t197) * t187, t1 * t65 - t12 * t206 - t152 * t2 + t207 * t7 - t208 * t8 + t3 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165 * t186, t164 * t186, -t96 ^ 2 - t97 ^ 2, t58 * t97 - t59 * t96 + t122, 0, 0, 0, 0, 0, t141, -t221, t141, -t212 - t220, t221, -t153 * t7 + t51 * t8 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, t212 - t220, t11, -t100 * t170 - t214 - t216, t169, -t153 * t66 + t150, -t123 * t9 + t51 * t66 + t149, -t17 * t51 + t150 + 0.2e1 * t163 - t211, pkin(4) * t21 - qJ(5) * t22 + (-t10 + t8) * t153 + (t7 - t200) * t51, 0.2e1 * t161 - t12 * t51 + t17 * t153 + (-0.2e1 * qJD(5) + t9) * t123 - t149, -pkin(4) * t2 + qJ(5) * t1 - t10 * t7 - t12 * t17 + t200 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169 + t210, t11, -t123 ^ 2 - t212, t123 * t8 + t2 + t211;];
tauc_reg = t13;
