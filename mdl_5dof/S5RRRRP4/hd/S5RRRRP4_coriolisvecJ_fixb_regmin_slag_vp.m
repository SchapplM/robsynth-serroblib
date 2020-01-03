% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP4
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
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:13
% EndTime: 2019-12-31 21:51:17
% DurationCPUTime: 1.19s
% Computational Cost: add. (2328->203), mult. (3941->260), div. (0->0), fcn. (2460->6), ass. (0->144)
t186 = pkin(8) + pkin(7);
t115 = cos(qJ(4));
t116 = cos(qJ(3));
t148 = qJD(4) * t115;
t151 = qJD(3) * t116;
t185 = -t115 * t151 - t116 * t148;
t109 = qJD(1) + qJD(2);
t112 = sin(qJ(4));
t113 = sin(qJ(3));
t77 = t112 * t116 + t115 * t113;
t68 = t77 * t109;
t108 = qJD(3) + qJD(4);
t117 = cos(qJ(2));
t165 = pkin(1) * qJD(1);
t143 = t117 * t165;
t89 = t186 * t113;
t106 = t116 * pkin(8);
t90 = t116 * pkin(7) + t106;
t57 = t112 * t90 + t115 * t89;
t156 = t115 * t116;
t158 = t112 * t113;
t76 = -t156 + t158;
t140 = qJD(3) * t186;
t78 = t113 * t140;
t79 = t116 * t140;
t172 = -t57 * qJD(4) - t112 * t79 - t115 * t78 + t76 * t143;
t184 = t172 * t108;
t183 = t68 ^ 2;
t114 = sin(qJ(2));
t98 = t114 * pkin(1) + pkin(7);
t181 = -pkin(8) - t98;
t141 = t109 * t158;
t66 = -t109 * t156 + t141;
t101 = -t116 * pkin(3) - pkin(2);
t70 = t101 * t109 - t143;
t22 = t66 * pkin(4) - t68 * qJ(5) + t70;
t125 = t108 * t158;
t50 = t125 + t185;
t166 = t185 * t109;
t33 = t109 * t125 + t166;
t51 = t108 * t77;
t34 = t51 * t109;
t152 = qJD(3) * t113;
t139 = t109 * t152;
t164 = pkin(1) * qJD(2);
t142 = qJD(1) * t164;
t95 = t114 * t142;
t71 = pkin(3) * t139 + t95;
t8 = t34 * pkin(4) + t33 * qJ(5) - t68 * qJD(5) + t71;
t180 = t22 * t50 - t8 * t77;
t179 = t22 * t51 + t8 * t76;
t178 = t117 * pkin(1);
t177 = t22 * t66;
t176 = t22 * t68;
t175 = t68 * t66;
t174 = t70 * t66;
t173 = t70 * t68;
t58 = -t112 * t89 + t115 * t90;
t171 = t58 * qJD(4) - t112 * t78 + t115 * t79 - t77 * t143;
t144 = t114 * t165;
t133 = t186 * t109 + t144;
t62 = t133 * t116;
t161 = t112 * t62;
t61 = t133 * t113;
t29 = -t115 * t61 - t161;
t170 = pkin(3) * t148 + qJD(5) - t29;
t169 = t70 * t51 + t71 * t76;
t168 = -t70 * t50 + t71 * t77;
t83 = -t109 * pkin(2) - t143;
t167 = t113 * t95 + t83 * t151;
t74 = t181 * t113;
t75 = t116 * t98 + t106;
t43 = t112 * t75 - t115 * t74;
t132 = qJD(3) * t181;
t145 = t117 * t164;
t59 = t113 * t132 + t116 * t145;
t60 = -t113 * t145 + t116 * t132;
t10 = -t43 * qJD(4) + t112 * t60 + t115 * t59;
t163 = t10 * t108;
t44 = t112 * t74 + t115 * t75;
t11 = t44 * qJD(4) + t112 * t59 - t115 * t60;
t162 = t11 * t108;
t160 = t115 * t62;
t159 = t109 * t113;
t157 = t114 * t116;
t118 = qJD(3) ^ 2;
t155 = t118 * t113;
t105 = t118 * t116;
t56 = qJD(3) * pkin(3) - t61;
t26 = t115 * t56 - t161;
t154 = qJD(5) - t26;
t153 = t113 ^ 2 - t116 ^ 2;
t150 = qJD(3) * t117;
t149 = qJD(4) * t112;
t147 = -qJD(1) - t109;
t146 = -qJD(2) + t109;
t103 = t114 * t164;
t102 = pkin(3) * t152;
t138 = t109 * t151;
t137 = t113 * t150;
t20 = -t108 * pkin(4) + t154;
t27 = t112 * t56 + t160;
t21 = t108 * qJ(5) + t27;
t104 = t108 * qJD(5);
t126 = qJD(3) * t133;
t129 = t117 * t142;
t41 = -t113 * t126 + t116 * t129;
t42 = -t113 * t129 - t116 * t126;
t130 = -t112 * t42 - t115 * t41 - t56 * t148 + t62 * t149;
t6 = t104 - t130;
t7 = t112 * t41 - t115 * t42 + t62 * t148 + t56 * t149;
t134 = -t20 * t50 - t21 * t51 - t6 * t76 + t7 * t77;
t131 = t171 * t108;
t17 = t51 * pkin(4) + t50 * qJ(5) - t77 * qJD(5) + t102;
t128 = -t17 + t144;
t28 = -t112 * t61 + t160;
t127 = pkin(3) * t149 - t28;
t35 = t68 * pkin(4) + t66 * qJ(5);
t124 = -t7 - t176;
t123 = t26 * t108 + t130;
t122 = t27 * t108 - t7;
t121 = -t144 + t102;
t120 = -t109 * t83 - t129;
t119 = -t114 * t159 + t116 * t150;
t45 = t76 * pkin(4) - t77 * qJ(5) + t101;
t107 = t109 ^ 2;
t100 = -pkin(2) - t178;
t99 = -t115 * pkin(3) - pkin(4);
t96 = t112 * pkin(3) + qJ(5);
t86 = t101 - t178;
t81 = 0.2e1 * t113 * t138;
t80 = t103 + t102;
t72 = t83 * t152;
t65 = -0.2e1 * t153 * t109 * qJD(3);
t47 = t51 * t108;
t46 = t50 * t108;
t40 = t45 - t178;
t30 = pkin(3) * t159 + t35;
t23 = -t66 ^ 2 + t183;
t18 = -t166 + (-t141 + t66) * t108;
t16 = t103 + t17;
t9 = -t33 * t77 - t68 * t50;
t1 = t33 * t76 - t77 * t34 + t50 * t66 - t68 * t51;
t2 = [0, 0, 0, 0, -t109 * t103 - t95, t147 * t145, t81, t65, t105, -t155, 0, t100 * t139 - t98 * t105 + t72 + (t147 * t157 - t137) * t164, t100 * t138 - t119 * t164 + t98 * t155 + t167, t9, t1, -t46, -t47, 0, t86 * t34 + t80 * t66 - t162 + t169, -t86 * t33 + t80 * t68 - t163 + t168, t16 * t66 + t40 * t34 - t162 + t179, -t10 * t66 + t11 * t68 - t43 * t33 - t44 * t34 + t134, -t16 * t68 + t40 * t33 + t163 + t180, t21 * t10 + t20 * t11 + t22 * t16 + t8 * t40 + t7 * t43 + t6 * t44; 0, 0, 0, 0, t109 * t144 - t95, t146 * t143, t81, t65, t105, -t155, 0, -pkin(2) * t139 - pkin(7) * t105 + t72 + (t146 * t157 + t137) * t165, -pkin(2) * t138 + pkin(7) * t155 + t119 * t165 + t167, t9, t1, -t46, -t47, 0, t101 * t34 + t121 * t66 - t131 + t169, -t101 * t33 + t121 * t68 + t168 - t184, -t128 * t66 + t45 * t34 - t131 + t179, t171 * t68 - t172 * t66 - t57 * t33 - t58 * t34 + t134, t128 * t68 + t45 * t33 + t180 + t184, -t128 * t22 + t171 * t20 + t172 * t21 + t8 * t45 + t7 * t57 + t6 * t58; 0, 0, 0, 0, 0, 0, -t113 * t107 * t116, t153 * t107, 0, 0, 0, t120 * t113, t120 * t116, t175, t23, t18, 0, 0, t28 * t108 - t173 + (-t108 * t149 - t66 * t159) * pkin(3) - t7, t29 * t108 + t174 + (-t108 * t148 - t68 * t159) * pkin(3) + t130, -t127 * t108 - t30 * t66 + t124, -t99 * t33 - t96 * t34 + (t127 + t21) * t68 + (t20 - t170) * t66, t170 * t108 + t30 * t68 - t177 + t6, t127 * t20 + t170 * t21 - t22 * t30 + t6 * t96 + t7 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, t23, t18, 0, 0, t122 - t173, t123 + t174, -t35 * t66 + t122 - t176, pkin(4) * t33 - t34 * qJ(5) + (t21 - t27) * t68 + (t20 - t154) * t66, t35 * t68 + 0.2e1 * t104 - t123 - t177, -t7 * pkin(4) + t6 * qJ(5) + t154 * t21 - t20 * t27 - t22 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, t18, -t108 ^ 2 - t183, -t21 * t108 - t124;];
tauc_reg = t2;
