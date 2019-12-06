% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:13
% EndTime: 2019-12-05 17:10:18
% DurationCPUTime: 1.53s
% Computational Cost: add. (2574->213), mult. (4955->300), div. (0->0), fcn. (3524->8), ass. (0->156)
t110 = sin(qJ(2));
t109 = sin(qJ(3));
t113 = cos(qJ(3));
t114 = cos(qJ(2));
t78 = t109 * t114 + t110 * t113;
t121 = t78 * qJD(2);
t160 = qJD(3) * t113;
t118 = (t110 * t160 + t121) * qJD(1);
t161 = qJD(3) * t109;
t162 = qJD(1) * t114;
t92 = qJD(2) * pkin(2) + t162;
t42 = t92 * t161 + t118;
t104 = qJD(2) + qJD(3);
t103 = qJD(4) + qJD(5);
t107 = sin(qJ(5));
t108 = sin(qJ(4));
t111 = cos(qJ(5));
t112 = cos(qJ(4));
t77 = t107 * t112 + t108 * t111;
t202 = t103 * t77;
t34 = t202 * t104;
t105 = t108 ^ 2;
t106 = t112 ^ 2;
t164 = t105 + t106;
t141 = qJD(2) * t162;
t163 = qJD(1) * t110;
t93 = t109 * t163;
t176 = t104 * t93;
t41 = (qJD(3) * t92 + t141) * t113 - t176;
t201 = t164 * t41;
t158 = qJD(4) * t112;
t167 = t111 * t112;
t199 = -qJD(5) * t167 - t111 * t158;
t168 = t107 * t108;
t75 = -t167 + t168;
t45 = t75 * t78;
t76 = t109 * t110 - t113 * t114;
t198 = t104 * t76;
t155 = pkin(2) * t160;
t72 = t113 * t162 - t93;
t197 = t155 - t72;
t71 = t78 * qJD(1);
t133 = pkin(2) * t161 - t71;
t66 = t109 * t92 + t113 * t163;
t64 = pkin(7) * t104 + t66;
t143 = pkin(8) * t104 + t64;
t128 = qJD(4) * t143;
t17 = -t108 * t128 + t112 * t41;
t46 = t143 * t108;
t43 = qJD(4) * pkin(4) - t46;
t196 = (qJD(5) * t43 + t17) * t111;
t115 = qJD(4) ^ 2;
t97 = pkin(2) * t109 + pkin(7);
t195 = t133 * t104 + t115 * t97;
t194 = -pkin(8) - pkin(7);
t193 = -pkin(8) - t97;
t192 = pkin(2) * t113;
t191 = pkin(3) * t104;
t190 = t41 * t78;
t189 = t42 * t76;
t65 = t113 * t92 - t93;
t99 = -pkin(4) * t112 - pkin(3);
t52 = t99 * t104 - t65;
t70 = t77 * t104;
t188 = t52 * t70;
t151 = t104 * t168;
t68 = -t104 * t167 + t151;
t187 = t70 * t68;
t73 = t193 * t108;
t101 = t112 * pkin(8);
t74 = t112 * t97 + t101;
t48 = -t107 * t74 + t111 * t73;
t142 = qJD(4) * t193;
t61 = t108 * t142 + t112 * t155;
t62 = -t108 * t155 + t112 * t142;
t186 = t48 * qJD(5) + t107 * t62 + t111 * t61 + t75 * t72;
t49 = t107 * t73 + t111 * t74;
t185 = -t49 * qJD(5) - t107 * t61 + t111 * t62 + t77 * t72;
t159 = qJD(4) * t108;
t147 = t104 * t159;
t32 = pkin(4) * t147 + t42;
t184 = t202 * t52 + t32 * t75;
t126 = t103 * t168;
t53 = t126 + t199;
t183 = t32 * t77 - t52 * t53;
t88 = t194 * t108;
t89 = pkin(7) * t112 + t101;
t59 = -t107 * t89 + t111 * t88;
t149 = qJD(4) * t194;
t79 = t108 * t149;
t80 = t112 * t149;
t182 = t59 * qJD(5) + t107 * t80 + t111 * t79 + t75 * t65;
t60 = t107 * t88 + t111 * t89;
t181 = -t60 * qJD(5) - t107 * t79 + t111 * t80 + t77 * t65;
t63 = -t65 - t191;
t179 = t42 * t108 + t63 * t158;
t154 = pkin(4) * t159;
t178 = t154 + t133;
t177 = t199 * t104;
t47 = t143 * t112;
t175 = t107 * t47;
t174 = t111 * t47;
t172 = t198 * t104;
t56 = t78 * qJD(3) + t121;
t171 = t56 * t104;
t170 = t66 * t104;
t169 = t104 * t108;
t166 = t115 * t108;
t165 = t105 - t106;
t157 = pkin(4) * t169;
t15 = t111 * t43 - t175;
t16 = t107 * t43 + t174;
t18 = -t108 * t41 - t112 * t128;
t138 = -qJD(5) * t175 + t107 * t18;
t3 = t138 + t196;
t139 = -t107 * t17 + t111 * t18;
t4 = -t16 * qJD(5) + t139;
t153 = t15 * t53 - t16 * t202 - t3 * t75 - t4 * t77;
t102 = t104 ^ 2;
t150 = t108 * t102 * t112;
t144 = -pkin(4) * t103 - t43;
t140 = -t104 * t63 - t41;
t137 = t164 * t65;
t134 = t112 * t147;
t132 = -t66 + t154;
t130 = pkin(7) * t115 - t170;
t129 = qJD(4) * (t65 - t191);
t127 = (-pkin(2) * t104 - t92) * qJD(3);
t125 = t115 * t78 + t171;
t124 = 0.2e1 * qJD(4) * t198;
t123 = t52 * t68 - t138;
t98 = -pkin(3) - t192;
t120 = qJD(4) * (t104 * t98 - t197);
t119 = t197 * t164;
t116 = qJD(2) ^ 2;
t100 = t115 * t112;
t86 = t99 - t192;
t83 = -0.2e1 * t134;
t82 = 0.2e1 * t134;
t67 = -0.2e1 * t165 * t104 * qJD(4);
t57 = t63 * t159;
t51 = t202 * t103;
t50 = t53 * t103;
t44 = t77 * t78;
t33 = t104 * t126 + t177;
t27 = -t68 ^ 2 + t70 ^ 2;
t22 = t70 * t103 - t34;
t21 = -t177 + (-t151 + t68) * t103;
t20 = -t111 * t46 - t175;
t19 = t107 * t46 - t174;
t11 = t202 * t68 + t34 * t75;
t10 = -t33 * t77 - t53 * t70;
t7 = t103 * t45 + t198 * t77;
t6 = t198 * t75 - t202 * t78;
t5 = -t202 * t70 + t33 * t75 - t34 * t77 + t53 * t68;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 * t110, -t116 * t114, 0, 0, 0, 0, 0, 0, 0, 0, -t171, t172, 0, -t198 * t66 - t56 * t65 + t189 + t190, 0, 0, 0, 0, 0, 0, t108 * t124 - t125 * t112, t125 * t108 + t112 * t124, -t164 * t172, t63 * t56 + t189 + t164 * (-t198 * t64 + t190), 0, 0, 0, 0, 0, 0, t103 * t7 + t34 * t76 + t56 * t68, -t103 * t6 - t33 * t76 + t56 * t70, -t33 * t44 + t34 * t45 - t6 * t68 - t7 * t70, t15 * t7 + t16 * t6 - t3 * t45 + t32 * t76 - t4 * t44 + t52 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t104 + t109 * t127 - t118, t72 * t104 + (t127 - t141) * t113 + t176, 0, t65 * t71 - t66 * t72 + (t109 * t41 - t113 * t42 + (-t109 * t65 + t113 * t66) * qJD(3)) * pkin(2), t82, t67, t100, t83, -t166, 0, t57 + t108 * t120 + (-t195 - t42) * t112, t195 * t108 + t112 * t120 + t179, t119 * t104 + t201, t119 * t64 + t133 * t63 + t201 * t97 + t42 * t98, t10, t5, -t50, t11, -t51, 0, t185 * t103 + t178 * t68 + t86 * t34 + t184, -t186 * t103 + t178 * t70 - t86 * t33 + t183, -t185 * t70 - t186 * t68 + t33 * t48 - t34 * t49 + t153, t185 * t15 + t186 * t16 + t178 * t52 + t3 * t49 + t32 * t86 + t4 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t170 - t42, t65 * t104 - t41, 0, 0, t82, t67, t100, t83, -t166, 0, t57 + t108 * t129 + (-t130 - t42) * t112, t130 * t108 + t112 * t129 + t179, -t104 * t137 + t201, -t42 * pkin(3) + pkin(7) * t201 - t64 * t137 - t63 * t66, t10, t5, -t50, t11, -t51, 0, t181 * t103 + t132 * t68 + t99 * t34 + t184, -t182 * t103 + t132 * t70 - t99 * t33 + t183, -t181 * t70 - t182 * t68 + t33 * t59 - t34 * t60 + t153, t132 * t52 + t181 * t15 + t182 * t16 + t3 * t60 + t32 * t99 + t4 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, t165 * t102, 0, t150, 0, 0, t140 * t108, t140 * t112, 0, 0, t187, t27, t21, -t187, t22, 0, -t68 * t157 - t103 * t19 - t188 + (t144 * t107 - t174) * qJD(5) + t139, -t70 * t157 + t103 * t20 + (t144 * qJD(5) - t17) * t111 + t123, (t16 + t19) * t70 + (-t15 + t20) * t68 + (-t107 * t34 + t111 * t33 + (t107 * t70 - t111 * t68) * qJD(5)) * pkin(4), -t15 * t19 - t16 * t20 + (-t52 * t169 + t107 * t3 + t111 * t4 + (-t107 * t15 + t111 * t16) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, t27, t21, -t187, t22, 0, t16 * t103 - t188 + t4, t15 * t103 + t123 - t196, 0, 0;];
tauc_reg = t1;
