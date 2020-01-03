% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:12
% EndTime: 2019-12-31 17:19:15
% DurationCPUTime: 0.99s
% Computational Cost: add. (1020->158), mult. (2266->257), div. (0->0), fcn. (1851->4), ass. (0->144)
t92 = sin(qJ(2));
t162 = qJ(4) * t92;
t94 = cos(qJ(2));
t110 = -t94 * pkin(2) - t92 * pkin(6);
t65 = -pkin(1) + t110;
t93 = cos(qJ(3));
t55 = t93 * t65;
t112 = -t93 * t162 + t55;
t91 = sin(qJ(3));
t174 = pkin(5) * t91;
t29 = (-pkin(3) - t174) * t94 + t112;
t170 = t94 * pkin(5);
t135 = t91 * t170;
t35 = t112 - t135;
t182 = -t29 + t35;
t178 = -t29 / 0.2e1;
t117 = t35 / 0.2e1 + t178;
t164 = -qJ(4) - pkin(6);
t67 = t164 * t93;
t181 = t117 * t67;
t138 = t92 * qJD(1);
t125 = t93 * t138;
t87 = t91 ^ 2;
t89 = t93 ^ 2;
t76 = t89 - t87;
t180 = t76 * qJD(2) - 0.2e1 * t91 * t125;
t161 = qJ(4) * t94;
t169 = t94 * pkin(6);
t172 = t92 * pkin(2);
t68 = -t169 + t172;
t61 = t93 * t68;
t167 = t91 * t92;
t82 = pkin(5) * t167;
t32 = t92 * pkin(3) - t93 * t161 + t61 + t82;
t177 = t32 / 0.2e1;
t175 = -t92 / 0.2e1;
t173 = t91 * pkin(3);
t171 = t94 * pkin(3);
t116 = pkin(5) + t173;
t62 = t116 * t92;
t168 = t62 * t91;
t166 = t92 * t93;
t88 = t92 ^ 2;
t165 = t93 * t88;
t75 = t89 + t87;
t90 = t94 ^ 2;
t77 = t90 - t88;
t163 = pkin(3) * qJD(3);
t126 = t171 / 0.2e1;
t3 = (t126 + t117) * t93;
t160 = t3 * qJD(1);
t133 = t93 * t170;
t42 = t91 * t65 + t133;
t36 = -t91 * t162 + t42;
t134 = pkin(5) * t166;
t60 = t91 * t68;
t37 = -t91 * t161 - t134 + t60;
t63 = t116 * t94;
t5 = t29 * t32 + t36 * t37 + t62 * t63;
t159 = t5 * qJD(1);
t6 = (t29 * t94 + t32 * t92) * t93 + (t36 * t94 + t37 * t92) * t91;
t158 = t6 * qJD(1);
t7 = t182 * t167;
t157 = t7 * qJD(1);
t136 = pkin(3) * t166;
t8 = t62 * t136 + t182 * t36;
t156 = t8 * qJD(1);
t155 = qJD(2) * t91;
t154 = qJD(2) * t92;
t153 = qJD(2) * t93;
t152 = qJD(2) * t94;
t151 = qJD(3) * t91;
t150 = qJD(3) * t93;
t149 = qJD(3) * t94;
t11 = (t29 * t93 + t36 * t91) * t92;
t148 = t11 * qJD(1);
t41 = -t55 + t135;
t12 = t41 * t92 + (-t82 + t61) * t94;
t147 = t12 * qJD(1);
t13 = t60 * t94 + (-t42 + t133) * t92;
t146 = t13 * qJD(1);
t30 = -t88 * t174 - t41 * t94;
t145 = t30 * qJD(1);
t31 = -pkin(5) * t165 - t42 * t94;
t144 = t31 * qJD(1);
t56 = t75 * t88;
t143 = t56 * qJD(1);
t58 = t77 * t91;
t142 = t58 * qJD(1);
t59 = t93 * t90 - t165;
t141 = t59 * qJD(1);
t140 = t75 * qJD(2);
t139 = t77 * qJD(1);
t137 = t94 * qJD(1);
t85 = -t93 * pkin(3) - pkin(2);
t132 = t85 * t166;
t131 = pkin(1) * t138;
t130 = pkin(1) * t137;
t129 = pkin(3) * t150;
t128 = pkin(3) * t151;
t127 = -t171 / 0.2e1;
t124 = t91 * t149;
t123 = t93 * t149;
t122 = t91 * t150;
t121 = t91 * t153;
t120 = t92 * t152;
t119 = t92 * t137;
t118 = t92 * t153;
t114 = t91 * t118;
t113 = -qJD(3) + t137;
t111 = -t67 * t175 + t178;
t1 = t181 + (t177 - t132 / 0.2e1 - t168 / 0.2e1) * pkin(3);
t14 = t85 * t173;
t109 = -t1 * qJD(1) + t14 * qJD(2);
t108 = t113 * t92;
t66 = t164 * t91;
t105 = (t36 / 0.2e1 + t66 * t175) * t93;
t10 = -t170 / 0.2e1 + t105 + (t127 + t111) * t91;
t39 = -t66 * t91 - t67 * t93;
t107 = t10 * qJD(1) + t39 * qJD(2);
t106 = t169 / 0.2e1 - t172 / 0.2e1;
t99 = t106 * t91;
t33 = t60 / 0.2e1 - t99;
t104 = pkin(2) * t153 - t33 * qJD(1);
t98 = t106 * t93;
t34 = -t61 / 0.2e1 + t98;
t103 = pkin(2) * t155 - t34 * qJD(1);
t102 = t93 * t108;
t50 = (t87 / 0.2e1 - t89 / 0.2e1) * t92;
t101 = -t50 * qJD(1) + t121;
t100 = t125 + t155;
t97 = t91 * qJD(1) * t165 + t50 * qJD(2);
t57 = t76 * t88;
t96 = t57 * qJD(1) + 0.2e1 * t114;
t83 = t154 / 0.2e1;
t54 = (t137 - qJD(3) / 0.2e1) * t92;
t49 = t100 * pkin(3);
t48 = t50 * qJD(3);
t23 = t82 + t61 / 0.2e1 + t98;
t22 = t134 - t60 / 0.2e1 - t99;
t9 = t170 / 0.2e1 + t105 + (t126 + t111) * t91;
t4 = (t117 + t127) * t93;
t2 = pkin(3) * t177 - t181 + (t132 + t168) * pkin(3) / 0.2e1;
t15 = [0, 0, 0, t120, t77 * qJD(2), 0, 0, 0, -pkin(1) * t154, -pkin(1) * t152, t120 * t89 - t122 * t88, -t57 * qJD(3) - 0.2e1 * t114 * t94, -t59 * qJD(2) + t124 * t92, t58 * qJD(2) + t123 * t92, -t120, -t12 * qJD(2) - t31 * qJD(3), t13 * qJD(2) + t30 * qJD(3), -t6 * qJD(2) - t7 * qJD(3) + t56 * qJD(4), t5 * qJD(2) + t8 * qJD(3) - t11 * qJD(4); 0, 0, 0, t119, t139, t152, -t154, 0, -pkin(5) * t152 - t131, pkin(5) * t154 - t130, -t48 + (t138 * t89 + t121) * t94, -0.2e1 * t92 * t122 + t180 * t94, t91 * t154 - t141, t118 + t142, -t54, -t147 + (t110 * t91 - t133) * qJD(2) + t23 * qJD(3), t146 + (t110 * t93 + t135) * qJD(2) + t22 * qJD(3), -t158 + ((-t66 * t94 + t37) * t93 + (t67 * t94 - t32) * t91) * qJD(2) + t4 * qJD(3), t159 + (t32 * t66 - t37 * t67 + t63 * t85) * qJD(2) + t2 * qJD(3) + t9 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96, t91 * t108, t102, t83, t23 * qJD(2) - t42 * qJD(3) - t144, t22 * qJD(2) + t41 * qJD(3) + t145, t4 * qJD(2) + t128 * t92 - t157, t2 * qJD(2) - t36 * t163 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t9 * qJD(2) - t148; 0, 0, 0, -t119, -t139, 0, 0, 0, t131, t130, -t119 * t89 - t48, 0.2e1 * t91 * t102, -t123 + t141, t124 - t142, t54, t34 * qJD(3) + t147, t33 * qJD(3) - t146, t3 * qJD(3) + t158, -t1 * qJD(3) + t10 * qJD(4) - t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t76 * qJD(3), 0, 0, 0, -pkin(2) * t151, -pkin(2) * t150, t75 * qJD(4), t14 * qJD(3) + t39 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t180, -t113 * t93, t113 * t91, -t138 / 0.2e1, -pkin(6) * t150 - t103, pkin(6) * t151 - t104, -t129 + t160, t67 * t163 + t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t96, (-t138 * t91 + t153) * t94, -t100 * t94, t83, -t34 * qJD(2) + t144, -t33 * qJD(2) - t145, -t3 * qJD(2) + t157, t1 * qJD(2) - qJD(4) * t136 - t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t180, t93 * t137, -t91 * t137, t138 / 0.2e1, t103, t104, -t160, -qJD(4) * t173 - t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, -t10 * qJD(2) + t129 * t92 + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, -t107 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t15;
