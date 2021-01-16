% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:46
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:45:23
% EndTime: 2021-01-15 16:45:31
% DurationCPUTime: 1.50s
% Computational Cost: add. (1580->264), mult. (4082->378), div. (0->0), fcn. (2849->8), ass. (0->147)
t87 = sin(qJ(3));
t139 = qJD(2) * t87;
t86 = sin(qJ(4));
t119 = t86 * t139;
t89 = cos(qJ(4));
t133 = t89 * qJD(3);
t60 = t119 - t133;
t90 = cos(qJ(3));
t131 = t90 * qJD(2);
t77 = -qJD(4) + t131;
t168 = t60 * t77;
t118 = t90 * t133;
t129 = qJD(3) * qJD(4);
t36 = -qJD(2) * t118 + qJD(4) * t119 - t89 * t129;
t184 = -t36 + t168;
t134 = t86 * qJD(3);
t62 = t89 * t139 + t134;
t167 = t62 * t77;
t135 = qJD(4) * t89;
t122 = t87 * t135;
t97 = t90 * t134 + t122;
t37 = t97 * qJD(2) + t86 * t129;
t183 = t37 - t167;
t84 = sin(pkin(5));
t142 = qJD(1) * t84;
t91 = cos(qJ(2));
t156 = t90 * t91;
t103 = pkin(3) * t87 - pkin(8) * t90;
t65 = t103 * qJD(3);
t88 = sin(qJ(2));
t182 = (-t86 * t156 + t88 * t89) * t142 - t87 * pkin(7) * t134 - t89 * t65;
t70 = -t90 * pkin(3) - t87 * pkin(8) - pkin(2);
t181 = -(t89 * t156 + t86 * t88) * t142 + t70 * t135 + t86 * t65;
t141 = qJD(1) * t90;
t121 = t88 * t142;
t68 = qJD(2) * pkin(7) + t121;
t85 = cos(pkin(5));
t41 = t85 * t141 - t87 * t68;
t180 = t62 ^ 2;
t179 = t60 * pkin(4);
t161 = t85 * t87;
t76 = qJD(1) * t161;
t42 = t90 * t68 + t76;
t34 = qJD(3) * pkin(8) + t42;
t120 = t91 * t142;
t44 = t70 * qJD(2) - t120;
t14 = -t86 * t34 + t89 * t44;
t10 = -t62 * qJ(5) + t14;
t5 = -t77 * pkin(4) + t10;
t178 = t10 - t5;
t144 = qJ(5) * t89;
t100 = pkin(4) * t87 - t90 * t144;
t132 = t89 * qJD(5);
t145 = qJ(5) * t87;
t157 = t89 * t90;
t78 = pkin(7) * t157;
t177 = t87 * t132 - t100 * qJD(3) - (-t78 + (-t70 + t145) * t86) * qJD(4) + t182;
t158 = t87 * t89;
t176 = -(-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t158 - (-qJD(5) * t87 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t90) * t86 - t181;
t160 = t86 * t44;
t15 = t89 * t34 + t160;
t11 = -t60 * qJ(5) + t15;
t175 = t11 * t77;
t105 = t87 * t120;
t137 = qJD(3) * t90;
t20 = qJD(2) * t105 + qJD(3) * t76 + t68 * t137;
t12 = t37 * pkin(4) + t20;
t174 = t12 * t86;
t173 = t12 * t89;
t172 = t20 * t86;
t171 = t20 * t89;
t33 = -qJD(3) * pkin(3) - t41;
t170 = t33 * t86;
t169 = t36 * t86;
t166 = t62 * t86;
t165 = t77 * t89;
t164 = t84 * t88;
t163 = t84 * t91;
t93 = qJD(2) ^ 2;
t162 = t84 * t93;
t159 = t86 * t90;
t92 = qJD(3) ^ 2;
t155 = t92 * t87;
t154 = t92 * t90;
t153 = -qJ(5) - pkin(8);
t111 = qJD(4) * t153;
t64 = t103 * qJD(2);
t113 = -t86 * t41 + t89 * t64;
t152 = t100 * qJD(2) + t86 * qJD(5) - t89 * t111 + t113;
t125 = t86 * t131;
t150 = t89 * t41 + t86 * t64;
t151 = -qJ(5) * t125 - t86 * t111 - t132 + t150;
t147 = t86 * t70 + t78;
t82 = t87 ^ 2;
t146 = -t90 ^ 2 + t82;
t143 = qJD(2) * pkin(2);
t140 = qJD(2) * t84;
t138 = qJD(3) * t87;
t136 = qJD(4) * t86;
t130 = qJD(2) * qJD(3);
t128 = t88 * t162;
t127 = t88 * t140;
t126 = t91 * t140;
t124 = t77 * t136;
t123 = t87 * t136;
t115 = t87 * t130;
t19 = -t68 * t138 + (qJD(3) * t85 + t126) * t141;
t40 = (t65 + t121) * qJD(2);
t114 = t86 * t19 - t89 * t40;
t112 = -qJD(5) - t179;
t110 = -t44 * t135 + t34 * t136 - t89 * t19 - t86 * t40;
t109 = t87 * t126;
t108 = t60 * t120;
t107 = t62 * t120;
t106 = t90 * t126;
t104 = pkin(4) * t115;
t69 = -t120 - t143;
t102 = -t69 - t120;
t101 = qJD(2) * t82 - t77 * t90;
t49 = t90 * t164 + t161;
t26 = -t89 * t163 - t49 * t86;
t99 = t86 * t163 - t49 * t89;
t48 = t87 * t164 - t85 * t90;
t98 = t37 * qJ(5) + t110;
t96 = qJD(3) * (-t102 - t143);
t95 = -t15 * qJD(4) - t114;
t94 = t36 * qJ(5) + t95;
t80 = -t89 * pkin(4) - pkin(3);
t72 = t153 * t89;
t71 = t153 * t86;
t66 = (pkin(4) * t86 + pkin(7)) * t87;
t59 = t89 * t70;
t57 = t60 ^ 2;
t43 = t97 * pkin(4) + pkin(7) * t137;
t28 = -t86 * t145 + t147;
t25 = t49 * qJD(3) + t109;
t24 = -t48 * qJD(3) + t106;
t22 = pkin(4) * t125 + t42;
t21 = -t87 * t144 + t59 + (-pkin(7) * t86 - pkin(4)) * t90;
t17 = -t112 + t33;
t8 = t26 * qJD(4) + t86 * t127 + t24 * t89;
t7 = t99 * qJD(4) + t89 * t127 - t24 * t86;
t4 = -t60 * qJD(5) - t98;
t3 = -t62 * qJD(5) + t104 + t94;
t2 = t115 * t99 + t25 * t62 - t48 * t36 + t8 * t77;
t1 = t26 * t115 + t25 * t60 + t48 * t37 - t7 * t77;
t6 = [0, 0, -t128, -t91 * t162, 0, 0, 0, 0, 0, -t90 * t128 + (-t25 - t109) * qJD(3), t87 * t128 + (-t24 - t106) * qJD(3), 0, 0, 0, 0, 0, t1, t2, t1, t2, t26 * t36 + t37 * t99 - t8 * t60 - t7 * t62, t11 * t8 + t12 * t48 + t17 * t25 + t3 * t26 - t4 * t99 + t5 * t7; 0, 0, 0, 0, 0.2e1 * t90 * t115, -0.2e1 * t146 * t130, t154, -t155, 0, -pkin(7) * t154 + t87 * t96, pkin(7) * t155 + t90 * t96, -t36 * t158 + (t118 - t123) * t62, (-t60 * t89 - t166) * t137 + (t169 - t37 * t89 + (t60 * t86 - t62 * t89) * qJD(4)) * t87, t77 * t123 + t36 * t90 + (t101 * t89 + t62 * t87) * qJD(3), t77 * t122 + t37 * t90 + (-t101 * t86 - t60 * t87) * qJD(3), (-t77 - t131) * t138, (t70 * t136 + t182) * t77 + ((pkin(7) * t60 + t170) * qJD(3) + (t160 + (pkin(7) * t77 + t34) * t89) * qJD(4) + t114) * t90 + (-t108 + t33 * t135 + pkin(7) * t37 + t172 + ((-pkin(7) * t159 + t59) * qJD(2) + t14) * qJD(3)) * t87, t181 * t77 + (t33 * t133 + (qJD(3) * t62 - t124) * pkin(7) - t110) * t90 + (-t107 - t33 * t136 - pkin(7) * t36 + t171 + (-pkin(7) * t165 - t147 * qJD(2) - t15) * qJD(3)) * t87, t66 * t37 + t43 * t60 + (t17 * t134 - t3) * t90 + t177 * t77 + (-t108 + t17 * t135 + t174 + (qJD(2) * t21 + t5) * qJD(3)) * t87, -t66 * t36 + t43 * t62 + (t17 * t133 + t4) * t90 - t176 * t77 + (-t107 - t17 * t136 + t173 + (-qJD(2) * t28 - t11) * qJD(3)) * t87, t21 * t36 - t28 * t37 + t177 * t62 + t176 * t60 + (-t11 * t86 - t5 * t89) * t137 + (-t3 * t89 - t4 * t86 + (-t11 * t89 + t5 * t86) * qJD(4)) * t87, t12 * t66 + t3 * t21 + t4 * t28 - t177 * t5 + (t43 - t105) * t17 - t176 * t11; 0, 0, 0, 0, -t87 * t93 * t90, t146 * t93, 0, 0, 0, t42 * qJD(3) - t69 * t139 - t20, t102 * t131, -t62 * t165 - t169, -t183 * t86 + t184 * t89, -t77 * t135 + (t77 * t157 + (-t62 + t134) * t87) * qJD(2), t124 + (-t77 * t159 + (t60 + t133) * t87) * qJD(2), t77 * t139, -pkin(3) * t37 - t171 + t113 * t77 - t42 * t60 + (pkin(8) * t165 + t170) * qJD(4) + (-t14 * t87 + (-pkin(8) * t138 - t33 * t90) * t86) * qJD(2), pkin(3) * t36 + t172 - t150 * t77 - t42 * t62 + (-t86 * pkin(8) * t77 + t33 * t89) * qJD(4) + (-t33 * t157 + (-pkin(8) * t133 + t15) * t87) * qJD(2), -t173 - t22 * t60 + t80 * t37 + t152 * t77 + (t17 + t179) * t136 + (-t17 * t159 + (qJD(3) * t71 - t5) * t87) * qJD(2), t174 - t22 * t62 - t80 * t36 - t151 * t77 + (pkin(4) * t166 + t17 * t89) * qJD(4) + (-t17 * t157 + (qJD(3) * t72 + t11) * t87) * qJD(2), t71 * t36 + t72 * t37 + t152 * t62 + t151 * t60 + (t77 * t5 + t4) * t89 + (-t3 + t175) * t86, t12 * t80 + t3 * t71 - t4 * t72 - t152 * t5 + (pkin(4) * t136 - t22) * t17 - t151 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60, -t57 + t180, -t36 - t168, -t167 - t37, t115, -t15 * t77 - t33 * t62 + t95, -t14 * t77 + t33 * t60 + t110, 0.2e1 * t104 - t175 + (t112 - t17) * t62 + t94, -t180 * pkin(4) - t10 * t77 + (qJD(5) + t17) * t60 + t98, t36 * pkin(4) + t178 * t60, -t178 * t11 + (-t17 * t62 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, t184, -t57 - t180, t11 * t60 + t5 * t62 + t12;];
tauc_reg = t6;
