% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:39
% EndTime: 2019-03-09 02:42:43
% DurationCPUTime: 1.50s
% Computational Cost: add. (1951->234), mult. (4695->313), div. (0->0), fcn. (3224->8), ass. (0->137)
t100 = sin(pkin(10));
t102 = cos(pkin(10));
t105 = sin(qJ(3));
t107 = cos(qJ(3));
t82 = t100 * t107 + t102 * t105;
t76 = t82 * qJD(1);
t177 = qJD(6) + t76;
t104 = sin(qJ(6));
t179 = t104 * t177;
t106 = cos(qJ(6));
t141 = qJD(1) * qJD(3);
t137 = t107 * t141;
t138 = t105 * t141;
t67 = -t100 * t138 + t102 * t137;
t63 = t106 * t67;
t121 = -t177 * t179 + t63;
t143 = t104 * qJD(3);
t146 = qJD(1) * t107;
t147 = qJD(1) * t105;
t73 = t100 * t147 - t102 * t146;
t54 = -t106 * t73 + t143;
t180 = t177 * t54;
t178 = qJD(6) - t177;
t92 = sin(pkin(9)) * pkin(1) + pkin(7);
t155 = qJ(4) + t92;
t71 = t76 ^ 2;
t176 = -t73 ^ 2 - t71;
t132 = t155 * qJD(1);
t57 = t107 * qJD(2) - t105 * t132;
t58 = t105 * qJD(2) + t107 * t132;
t49 = t100 * t58;
t27 = t102 * t57 - t49;
t150 = -qJD(5) + t27;
t169 = t73 * pkin(5);
t159 = t102 * t58;
t26 = t100 * t57 + t159;
t93 = -pkin(3) * t102 - pkin(4);
t88 = -pkin(8) + t93;
t52 = qJD(3) * pkin(3) + t57;
t24 = t100 * t52 + t159;
t21 = -qJD(3) * qJ(5) - t24;
t9 = -t21 - t169;
t174 = t88 * t67 + (t9 - t26 + t169) * t177;
t144 = qJD(6) * t106;
t153 = t102 * t107;
t81 = t100 * t105 - t153;
t139 = t81 * t144;
t162 = t81 * t67;
t75 = t82 * qJD(3);
t173 = -t104 * (t177 * t75 + t162) - t177 * t139;
t171 = pkin(4) + pkin(8);
t66 = qJD(1) * t75;
t170 = t66 * pkin(4);
t168 = t76 * pkin(5);
t140 = qJD(1) * qJD(4);
t46 = qJD(3) * t57 + t107 * t140;
t47 = -qJD(3) * t58 - t105 * t140;
t10 = t100 * t46 - t102 * t47;
t79 = t155 * t105;
t80 = t155 * t107;
t42 = t100 * t80 + t102 * t79;
t167 = t10 * t42;
t166 = t10 * t81;
t94 = -cos(pkin(9)) * pkin(1) - pkin(2);
t123 = -t107 * pkin(3) + t94;
t113 = -t82 * qJ(5) + t123;
t28 = t171 * t81 + t113;
t165 = t28 * t67;
t56 = qJD(3) * t106 + t104 * t73;
t164 = t56 * t73;
t163 = t73 * t54;
t34 = -qJD(6) * t143 + t104 * t66 + t144 * t73;
t145 = qJD(3) * t105;
t78 = qJD(3) * t153 - t100 * t145;
t161 = t34 * t82 + t56 * t78;
t11 = t100 * t47 + t102 * t46;
t160 = t105 ^ 2 - t107 ^ 2;
t157 = t106 * t177;
t156 = t34 * t106;
t84 = qJD(1) * t94;
t154 = qJD(6) * t81;
t108 = qJD(3) ^ 2;
t152 = t108 * t105;
t151 = t108 * t107;
t149 = t168 - t150;
t96 = pkin(3) * t145;
t89 = pkin(3) * t138;
t136 = -t67 * qJ(5) + t89;
t135 = pkin(3) * t147 + t73 * qJ(5);
t134 = qJD(3) * t155;
t59 = t107 * qJD(4) - t105 * t134;
t60 = -t105 * qJD(4) - t107 * t134;
t29 = t100 * t59 - t102 * t60;
t23 = t102 * t52 - t49;
t130 = qJD(5) - t23;
t62 = t106 * t66;
t35 = qJD(6) * t56 - t62;
t128 = -t35 * t82 - t54 * t78;
t8 = -qJD(3) * qJD(5) - t11;
t115 = t123 * qJD(1);
t72 = qJD(4) + t115;
t111 = -t76 * qJ(5) + t72;
t17 = t171 * t73 + t111;
t7 = -qJD(3) * t171 + t130 + t168;
t2 = t104 * t7 + t106 * t17;
t126 = t104 * t17 - t106 * t7;
t30 = t100 * t60 + t102 * t59;
t43 = -t100 * t79 + t102 * t80;
t125 = 0.2e1 * qJD(3) * t84;
t31 = t73 * pkin(4) + t111;
t122 = t31 * t76 + t10;
t120 = -t76 * qJD(5) + t136;
t119 = -t78 * qJ(5) - t82 * qJD(5) + t96;
t3 = -pkin(5) * t66 - t8;
t118 = t3 + (-qJD(6) * t88 + t171 * t76 + t135) * t177;
t117 = -t154 * t179 + t157 * t75 + t63 * t81;
t32 = pkin(5) * t82 + t42;
t116 = t3 * t81 - t32 * t67 + t9 * t75;
t114 = -t104 * t67 - t157 * t177;
t112 = -t66 * t82 - t73 * t78 + t75 * t76 + t162;
t110 = t10 * t82 + t29 * t76 - t30 * t73 + t42 * t67 - t43 * t66;
t109 = qJD(1) ^ 2;
t90 = pkin(3) * t100 + qJ(5);
t69 = qJD(3) * t73;
t38 = pkin(4) * t81 + t113;
t36 = pkin(4) * t76 + t135;
t33 = -pkin(5) * t81 + t43;
t25 = pkin(4) * t75 + t119;
t19 = -qJD(3) * pkin(4) + t130;
t18 = t120 + t170;
t16 = -pkin(5) * t75 + t30;
t15 = pkin(5) * t78 + t29;
t12 = t171 * t75 + t119;
t6 = t171 * t66 + t120;
t5 = pkin(5) * t67 + t10;
t4 = t106 * t5;
t1 = [0, 0, 0, 0, 0.2e1 * t105 * t137, -0.2e1 * t160 * t141, t151, -t152, 0, t105 * t125 - t151 * t92, t107 * t125 + t152 * t92, -t11 * t81 - t23 * t78 - t24 * t75 + t110, t167 + t11 * t43 - t23 * t29 + t24 * t30 + (t72 + t115) * t96, t19 * t78 + t21 * t75 + t8 * t81 + t110, qJD(3) * t29 - t18 * t81 - t25 * t73 - t31 * t75 - t38 * t66, qJD(3) * t30 - t18 * t82 - t25 * t76 - t31 * t78 - t38 * t67, t18 * t38 + t19 * t29 - t21 * t30 + t25 * t31 - t43 * t8 + t167, t56 * t139 + (t34 * t81 + t56 * t75) * t104 (-t104 * t54 + t106 * t56) * t75 + (-t104 * t35 + t156 + (-t104 * t56 - t106 * t54) * qJD(6)) * t81, t161 - t173, t117 + t128, t177 * t78 + t67 * t82, -t126 * t78 + t16 * t54 + t33 * t35 + t4 * t82 + (-t12 * t177 - t6 * t82 - t165) * t104 + (t15 * t177 - t116) * t106 + ((-t104 * t32 - t106 * t28) * t177 - t2 * t82 + t9 * t104 * t81) * qJD(6), t16 * t56 - t2 * t78 + t33 * t34 + (-(qJD(6) * t32 + t12) * t177 - t165 - (qJD(6) * t7 + t6) * t82 + t9 * t154) * t106 + (-(-qJD(6) * t28 + t15) * t177 - (-qJD(6) * t17 + t5) * t82 + t116) * t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, -t151, t112, t11 * t82 - t23 * t75 + t24 * t78 + t166, t112, t75 * qJD(3), t78 * qJD(3), t19 * t75 - t21 * t78 - t8 * t82 + t166, 0, 0, 0, 0, 0, t117 - t128, t161 + t173; 0, 0, 0, 0, -t105 * t109 * t107, t160 * t109, 0, 0, 0, -t84 * t147, -t84 * t146 (t24 - t26) * t76 + (-t23 + t27) * t73 + (-t100 * t66 - t102 * t67) * pkin(3), t23 * t26 - t24 * t27 + (-t10 * t102 + t100 * t11 - t147 * t72) * pkin(3), -t90 * t66 + t93 * t67 + (-t21 - t26) * t76 + (t19 + t150) * t73, -qJD(3) * t26 + t36 * t73 + t122, -t31 * t73 + t36 * t76 + (0.2e1 * qJD(5) - t27) * qJD(3) + t11, t10 * t93 + t150 * t21 - t19 * t26 - t31 * t36 - t8 * t90, -t179 * t56 + t156 (-t177 * t56 - t35) * t106 + (-t34 + t180) * t104, t121 + t164, t114 - t163, t177 * t73, t104 * t118 + t106 * t174 - t126 * t73 + t149 * t54 + t90 * t35, -t104 * t174 + t106 * t118 + t149 * t56 - t2 * t73 + t90 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t23 * t76 + t24 * t73 + t89, t176, -0.2e1 * t76 * qJD(3), -t67 + t69, t170 - t21 * t73 + (-qJD(5) - t19) * t76 + t136, 0, 0, 0, 0, 0, t114 + t163, -t121 + t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 + t69, -t76 * t73, -t71 - t108, qJD(3) * t21 + t122, 0, 0, 0, 0, 0, -qJD(3) * t54 + t121, -qJD(3) * t56 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t54 ^ 2 + t56 ^ 2, t34 + t180, -t178 * t56 + t62, t67, -t104 * t6 - t178 * t2 - t9 * t56 + t4, -t104 * t5 - t106 * t6 + t126 * t178 + t9 * t54;];
tauc_reg  = t1;
