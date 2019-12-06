% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRR5
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:28
% EndTime: 2019-12-05 18:16:33
% DurationCPUTime: 1.05s
% Computational Cost: add. (2559->192), mult. (5056->252), div. (0->0), fcn. (3111->8), ass. (0->143)
t112 = sin(qJ(3));
t115 = cos(qJ(3));
t182 = sin(pkin(9)) * pkin(1);
t147 = qJD(1) * t182;
t96 = cos(pkin(9)) * pkin(1) + pkin(2);
t86 = t96 * qJD(1);
t67 = -t112 * t147 + t115 * t86;
t113 = cos(qJ(5));
t114 = cos(qJ(4));
t150 = qJD(4) * t114;
t155 = t113 * t114;
t188 = -qJD(5) * t155 - t113 * t150;
t105 = qJD(1) + qJD(3);
t110 = sin(qJ(5));
t111 = sin(qJ(4));
t78 = t110 * t114 + t113 * t111;
t74 = t78 * t105;
t151 = qJD(4) * t111;
t100 = t114 * qJD(2);
t64 = t67 * qJD(3);
t171 = qJD(4) * t100 + t114 * t64;
t68 = t112 * t86 + t115 * t147;
t58 = t105 * pkin(7) + t68;
t22 = -t58 * t151 + t171;
t165 = t111 * t64;
t149 = t111 * qJD(2);
t46 = t114 * t58 + t149;
t23 = -t46 * qJD(4) - t165;
t166 = t111 * t58;
t45 = t100 - t166;
t187 = -t23 * t111 + t22 * t114 + (-t111 * t46 - t114 * t45) * qJD(4);
t133 = -t112 * t182 + t115 * t96;
t137 = pkin(8) * t105 + t58;
t128 = t137 * t111;
t17 = -qJD(4) * t128 + t171;
t43 = t100 - t128;
t40 = qJD(4) * pkin(4) + t43;
t186 = (qJD(5) * t40 + t17) * t113;
t185 = t68 * qJD(3);
t104 = qJD(4) + qJD(5);
t44 = t137 * t114 + t149;
t184 = -pkin(8) - pkin(7);
t169 = t112 * t96 + t115 * t182;
t76 = pkin(7) + t169;
t183 = -pkin(8) - t76;
t181 = t105 * pkin(3);
t180 = t114 * pkin(4);
t99 = -pkin(3) - t180;
t47 = t99 * t105 - t67;
t179 = t47 * t74;
t156 = t110 * t111;
t144 = t105 * t156;
t72 = -t105 * t155 + t144;
t178 = t74 * t72;
t141 = t105 * t151;
t50 = pkin(4) * t141 + t185;
t54 = t104 * t78;
t77 = -t155 + t156;
t177 = t47 * t54 + t50 * t77;
t89 = t184 * t111;
t102 = t114 * pkin(8);
t90 = t114 * pkin(7) + t102;
t62 = -t110 * t90 + t113 * t89;
t142 = qJD(4) * t184;
t79 = t111 * t142;
t80 = t114 * t142;
t176 = t62 * qJD(5) + t110 * t80 + t113 * t79 + t77 * t67;
t63 = t110 * t89 + t113 * t90;
t175 = -t63 * qJD(5) - t110 * t79 + t113 * t80 + t78 * t67;
t42 = t54 * t105;
t125 = t104 * t156;
t53 = t125 + t188;
t174 = -t78 * t42 + t53 * t72;
t173 = -t47 * t53 + t50 * t78;
t57 = -t67 - t181;
t172 = t111 * t185 + t57 * t150;
t170 = t188 * t105;
t168 = t105 * t57;
t167 = t110 * t44;
t164 = t113 * t44;
t48 = t53 * t104;
t161 = t67 * t105;
t160 = t68 * t105;
t69 = t133 * qJD(3);
t159 = t69 * t105;
t70 = t169 * qJD(3);
t158 = t70 * t105;
t157 = t105 * t111;
t116 = qJD(4) ^ 2;
t154 = t116 * t111;
t106 = t111 ^ 2;
t107 = t114 ^ 2;
t153 = t106 - t107;
t152 = t106 + t107;
t148 = pkin(4) * t157;
t146 = pkin(4) * t151;
t12 = t113 * t40 - t167;
t13 = t110 * t40 + t164;
t18 = -t44 * qJD(4) - t165;
t134 = -qJD(5) * t167 + t110 * t18;
t3 = t134 + t186;
t135 = -t110 * t17 + t113 * t18;
t4 = -t13 * qJD(5) + t135;
t145 = t12 * t53 - t13 * t54 - t3 * t77 - t4 * t78;
t103 = t105 ^ 2;
t143 = t111 * t103 * t114;
t138 = -pkin(4) * t104 - t40;
t136 = qJD(4) * t183;
t75 = -pkin(3) - t133;
t131 = t114 * t141;
t130 = -t68 + t146;
t41 = t105 * t125 + t170;
t129 = -t77 * t41 + t74 * t54;
t127 = pkin(7) * t116 - t160;
t126 = qJD(4) * (t67 - t181);
t124 = t116 * t76 + t158;
t59 = t183 * t111;
t60 = t114 * t76 + t102;
t27 = -t110 * t60 + t113 * t59;
t28 = t110 * t59 + t113 * t60;
t123 = t111 * t45 - t114 * t46;
t122 = qJD(4) * (t105 * t75 - t69);
t121 = t47 * t72 - t134;
t101 = t116 * t114;
t82 = -0.2e1 * t131;
t81 = 0.2e1 * t131;
t71 = -0.2e1 * t153 * t105 * qJD(4);
t66 = t75 - t180;
t61 = t70 + t146;
t51 = t57 * t151;
t49 = t54 * t104;
t35 = -t111 * t69 + t114 * t136;
t34 = t111 * t136 + t114 * t69;
t26 = -t72 ^ 2 + t74 ^ 2;
t20 = -t170 + (-t144 + t72) * t104;
t15 = t113 * t43 - t167;
t14 = -t110 * t43 - t164;
t11 = t42 * t77 + t72 * t54;
t10 = -t41 * t78 - t74 * t53;
t7 = -t28 * qJD(5) - t110 * t34 + t113 * t35;
t6 = t27 * qJD(5) + t110 * t35 + t113 * t34;
t5 = -t129 + t174;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185 - t158, -t64 - t159, 0, -t133 * t185 + t64 * t169 - t67 * t70 + t68 * t69, t81, t71, t101, t82, -t154, 0, t51 + t111 * t122 + (-t124 - t185) * t114, t124 * t111 + t114 * t122 + t172, t152 * t159 + t187, -t123 * t69 + t185 * t75 + t187 * t76 + t57 * t70, t10, t5, -t48, t11, -t49, 0, t7 * t104 + t66 * t42 + t61 * t72 + t177, -t6 * t104 - t66 * t41 + t61 * t74 + t173, t27 * t41 - t28 * t42 - t6 * t72 - t7 * t74 + t145, t12 * t7 + t13 * t6 + t4 * t27 + t3 * t28 + t47 * t61 + t50 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, -t101, 0, -t123 * qJD(4) + t22 * t111 + t23 * t114, 0, 0, 0, 0, 0, 0, -t49, t48, t129 + t174, -t12 * t54 - t13 * t53 + t3 * t78 - t4 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185 + t160, -t64 + t161, 0, 0, t81, t71, t101, t82, -t154, 0, t51 + t111 * t126 + (-t127 - t185) * t114, t127 * t111 + t114 * t126 + t172, -t152 * t161 + t187, -pkin(3) * t185 + pkin(7) * t187 + t123 * t67 - t57 * t68, t10, t5, -t48, t11, -t49, 0, t175 * t104 + t130 * t72 + t99 * t42 + t177, -t176 * t104 + t130 * t74 - t99 * t41 + t173, -t175 * t74 - t176 * t72 + t62 * t41 - t63 * t42 + t145, t175 * t12 + t176 * t13 + t130 * t47 + t3 * t63 + t4 * t62 + t50 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t153 * t103, 0, t143, 0, 0, (-t64 - t168) * t111, -t114 * t168 + (t45 + t166) * qJD(4) - t171, 0, 0, t178, t26, t20, -t178, 0, 0, -t72 * t148 - t14 * t104 - t179 + (t138 * t110 - t164) * qJD(5) + t135, -t74 * t148 + t15 * t104 + (t138 * qJD(5) - t17) * t113 + t121, (t13 + t14) * t74 + (-t12 + t15) * t72 + (-t110 * t42 + t113 * t41 + (t110 * t74 - t113 * t72) * qJD(5)) * pkin(4), -t12 * t14 - t13 * t15 + (-t47 * t157 + t110 * t3 + t113 * t4 + (-t110 * t12 + t113 * t13) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, t26, t20, -t178, 0, 0, t13 * t104 - t179 + t4, t12 * t104 + t121 - t186, 0, 0;];
tauc_reg = t1;
