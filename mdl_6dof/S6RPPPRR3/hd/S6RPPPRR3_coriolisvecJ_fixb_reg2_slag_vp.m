% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPPRR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:59
% EndTime: 2019-03-09 01:34:04
% DurationCPUTime: 2.03s
% Computational Cost: add. (4293->265), mult. (8578->370), div. (0->0), fcn. (5775->8), ass. (0->147)
t185 = cos(qJ(5));
t99 = cos(pkin(10));
t141 = t185 * t99;
t103 = sin(qJ(5));
t97 = sin(pkin(10));
t162 = t103 * t97;
t110 = t141 - t162;
t100 = cos(pkin(9));
t148 = qJD(1) * t100;
t72 = t103 * t99 + t185 * t97;
t70 = t72 * qJD(5);
t98 = sin(pkin(9));
t169 = t110 * t148 + t98 * t70;
t102 = sin(qJ(6));
t104 = cos(qJ(6));
t137 = qJD(1) * t162;
t79 = qJD(1) * t141;
t67 = -t79 + t137;
t149 = qJD(6) - t67;
t132 = t104 * t149;
t61 = qJD(1) * t70;
t195 = -t102 * t61 + t132 * t149;
t77 = qJD(5) * t79;
t60 = qJD(5) * t137 - t77;
t143 = qJD(1) * qJD(2);
t82 = t98 * t143;
t28 = -t61 * pkin(5) - t60 * pkin(8) + t82;
t105 = -pkin(1) - pkin(2);
t78 = t105 * qJD(1) + qJD(2);
t63 = qJ(2) * t148 + t98 * t78;
t51 = -qJD(1) * qJ(4) + t63;
t90 = t99 * qJD(3);
t31 = t90 + (pkin(7) * qJD(1) - t51) * t97;
t156 = qJD(1) * t99;
t38 = t97 * qJD(3) + t99 * t51;
t32 = -pkin(7) * t156 + t38;
t16 = t103 * t31 + t185 * t32;
t14 = qJD(5) * pkin(8) + t16;
t151 = t98 * qJD(1);
t62 = -qJ(2) * t151 + t100 * t78;
t48 = qJD(1) * pkin(3) + qJD(4) - t62;
t42 = pkin(4) * t156 + t48;
t68 = qJD(1) * t72;
t19 = -t67 * pkin(5) + t68 * pkin(8) + t42;
t6 = t102 * t19 + t104 * t14;
t112 = t103 * t32 - t185 * t31;
t145 = t100 * qJD(2);
t81 = -qJD(4) + t145;
t75 = t81 * qJD(1);
t193 = t110 * t75;
t9 = -t112 * qJD(5) + t193;
t2 = -qJD(6) * t6 - t102 * t9 + t104 * t28;
t190 = t149 * t6 + t2;
t115 = t102 * t14 - t104 * t19;
t1 = -t115 * qJD(6) + t102 * t28 + t104 * t9;
t127 = t115 * t149 + t1;
t194 = t72 * t75;
t167 = t97 ^ 2 + t99 ^ 2;
t136 = t167 * t75;
t191 = t110 * qJD(5);
t47 = t102 * qJD(5) - t104 * t68;
t155 = qJD(6) * t47;
t27 = t102 * t60 + t155;
t174 = t72 * t61;
t119 = -t149 * t191 + t174;
t146 = qJD(6) * t104;
t138 = t72 * t146;
t189 = -t102 * t119 + t138 * t149;
t188 = t68 ^ 2;
t187 = 0.2e1 * t82;
t168 = t100 * qJ(2) + t98 * t105;
t73 = -qJ(4) + t168;
t186 = pkin(7) - t73;
t10 = t16 * qJD(5) + t194;
t56 = t186 * t97;
t57 = t186 * t99;
t111 = t103 * t57 + t185 * t56;
t184 = t10 * t111;
t58 = t72 * t98;
t183 = t10 * t58;
t182 = t10 * t110;
t181 = t10 * t72;
t144 = t104 * qJD(5);
t45 = -t102 * t68 - t144;
t180 = t45 * t67;
t179 = t47 * t45;
t178 = t47 * t68;
t177 = t61 * t110;
t176 = t67 * t68;
t175 = t68 * t45;
t59 = t110 * t98;
t39 = -t104 * t100 - t102 * t59;
t173 = t39 * qJD(6) - t102 * t151 - t169 * t104;
t114 = t102 * t100 - t104 * t59;
t172 = t114 * qJD(6) + t169 * t102 - t104 * t151;
t171 = -t102 * t27 - t45 * t146;
t170 = t100 * t68 - t191 * t98;
t166 = t102 * t45;
t165 = t102 * t47;
t161 = t104 * t45;
t160 = t104 * t47;
t147 = qJD(6) * t102;
t26 = -qJD(6) * t144 - t104 * t60 - t68 * t147;
t159 = t26 * t102;
t158 = t27 * t104;
t106 = qJD(1) ^ 2;
t157 = t98 * t106;
t154 = t100 * t106;
t153 = t191 * qJD(5);
t152 = t70 * qJD(5);
t150 = t98 * qJD(2);
t142 = 0.2e1 * t143;
t139 = t72 * t147;
t135 = -0.2e1 * t68;
t134 = -t98 * qJ(2) + t100 * t105;
t133 = t102 * t149;
t129 = -pkin(8) * qJD(6) * t149 - t10;
t128 = pkin(3) - t134;
t13 = -qJD(5) * pkin(5) + t112;
t126 = -t13 * t191 - t181;
t125 = -t102 * t6 + t104 * t115;
t124 = -t102 * t115 - t104 * t6;
t123 = t110 * t26 + t70 * t47;
t122 = -t110 * t27 + t70 * t45;
t121 = -(-t97 * t51 + t90) * t97 + t38 * t99;
t120 = -t110 * t60 - t68 * t70;
t118 = -t191 * t67 - t174;
t117 = -t63 * t100 + t62 * t98;
t25 = t103 * t56 - t185 * t57;
t65 = t99 * pkin(4) + t128;
t29 = pkin(5) * t110 + t72 * pkin(8) + t65;
t12 = t102 * t29 + t104 * t25;
t11 = -t102 * t25 + t104 * t29;
t113 = -t104 * t61 + (t102 * t67 - t147) * t149;
t109 = pkin(8) * t61 + t149 * t13;
t108 = t119 * t104 + t139 * t149;
t107 = t125 * qJD(6) + t1 * t104 - t2 * t102;
t66 = t67 ^ 2;
t34 = -t68 * pkin(5) - t67 * pkin(8);
t30 = -t70 * pkin(5) + pkin(8) * t191 + t150;
t18 = t25 * qJD(5) + t72 * t81;
t17 = t111 * qJD(5) + t110 * t81;
t8 = t102 * t34 - t104 * t112;
t7 = t102 * t112 + t104 * t34;
t4 = -t12 * qJD(6) - t102 * t17 + t104 * t30;
t3 = t11 * qJD(6) + t102 * t30 + t104 * t17;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, qJ(2) * t142, 0, 0, 0, 0, 0, 0, t187, t100 * t142, 0 ((t100 * t168 - t98 * t134) * qJD(1) - t117) * qJD(2), 0, 0, 0, 0, 0, 0, t99 * t187, -0.2e1 * t97 * t82, -0.2e1 * t136, t121 * t81 + t73 * t136 + (qJD(1) * t128 + t48) * t150, t191 * t68 - t60 * t72, t118 + t120, -t153, t67 * t70 - t177, t152, 0, -t18 * qJD(5) - t42 * t70 - t65 * t61 + (qJD(1) * t110 - t67) * t150, -t17 * qJD(5) + t135 * t150 - t191 * t42 + t65 * t60, -t110 * t9 - t111 * t60 - t112 * t191 + t16 * t70 + t17 * t67 - t18 * t68 + t25 * t61 - t181, -t184 + t112 * t18 + t16 * t17 + t9 * t25 + (qJD(1) * t65 + t42) * t150, t47 * t139 + (-t191 * t47 + t26 * t72) * t104 -(-t161 - t165) * t191 + (-t159 + t158 + (t160 - t166) * qJD(6)) * t72, t108 - t123, -t45 * t138 + (-t191 * t45 - t27 * t72) * t102, t122 + t189, -t149 * t70 - t177, t102 * t126 - t11 * t61 + t110 * t2 - t111 * t27 + t115 * t70 - t13 * t138 + t149 * t4 + t18 * t45, -t1 * t110 + t104 * t126 + t111 * t26 + t12 * t61 + t13 * t139 - t149 * t3 + t18 * t47 + t6 * t70, t11 * t26 - t12 * t27 - t3 * t45 - t4 * t47 - t125 * t191 + (-qJD(6) * t124 + t1 * t102 + t2 * t104) * t72, t1 * t12 + t2 * t11 - t115 * t4 + t13 * t18 + t6 * t3 - t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t106 * qJ(2), 0, 0, 0, 0, 0, 0, -t157, -t154, 0, t117 * qJD(1), 0, 0, 0, 0, 0, 0, -t99 * t157, t97 * t157, t167 * t154, t98 * t136 + (-t48 * t98 + (-t121 - t150) * t100) * qJD(1), 0, 0, 0, 0, 0, 0, t170 * qJD(5) + t100 * t61 + t67 * t151, t169 * qJD(5) - t100 * t60 + t68 * t151, -t169 * t67 + t170 * t68 + t58 * t60 + t59 * t61, t183 + t9 * t59 - t169 * t16 - t170 * t112 + (-t42 - t145) * t151, 0, 0, 0, 0, 0, 0, t149 * t172 - t170 * t45 + t58 * t27 - t39 * t61, -t114 * t61 - t149 * t173 - t170 * t47 - t58 * t26, t114 * t27 - t172 * t47 - t173 * t45 + t39 * t26, -t1 * t114 - t115 * t172 - t170 * t13 + t173 * t6 + t2 * t39 + t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, -t153, -t118 + t120, t112 * t70 + t16 * t191 + t9 * t72 - t182, 0, 0, 0, 0, 0, 0, t122 - t189, t108 + t123 -(t161 - t165) * t191 + (-t159 - t158 + (t160 + t166) * qJD(6)) * t72, t107 * t72 - t124 * t191 + t13 * t70 - t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167 * t106, qJD(1) * t121 + t82, 0, 0, 0, 0, 0, 0, t135 * qJD(5), -t77 + (t67 + t137) * qJD(5), -t66 - t188, t112 * t68 - t16 * t67 + t82, 0, 0, 0, 0, 0, 0, t113 + t175, t178 - t195 (t26 + t180) * t104 + t47 * t133 + t171, t127 * t102 + t190 * t104 + t13 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, -t66 + t188, -t77 + (-t67 + t137) * qJD(5), -t176, 0, 0, t42 * t68 - t194, -t42 * t67 - t193, 0, 0, t132 * t47 - t159 (-t26 + t180) * t104 - t149 * t165 + t171, t178 + t195, t133 * t45 - t158, t113 - t175, t149 * t68, -pkin(5) * t27 + t102 * t109 + t104 * t129 - t115 * t68 - t149 * t7 - t16 * t45, pkin(5) * t26 - t102 * t129 + t104 * t109 + t149 * t8 - t16 * t47 - t6 * t68, t8 * t45 + t7 * t47 + ((-t27 + t155) * pkin(8) + t127) * t104 + ((qJD(6) * t45 - t26) * pkin(8) - t190) * t102, -t10 * pkin(5) + pkin(8) * t107 + t115 * t7 - t13 * t16 - t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t45 ^ 2 + t47 ^ 2, t149 * t45 - t26, -t179, t149 * t47 - t27, -t61, -t13 * t47 + t190, t13 * t45 - t127, 0, 0;];
tauc_reg  = t5;
