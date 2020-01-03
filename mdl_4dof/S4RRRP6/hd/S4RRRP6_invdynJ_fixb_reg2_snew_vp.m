% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:19:10
% EndTime: 2019-12-31 17:19:14
% DurationCPUTime: 1.02s
% Computational Cost: add. (2410->194), mult. (4918->236), div. (0->0), fcn. (3093->6), ass. (0->135)
t117 = sin(qJ(3));
t120 = cos(qJ(3));
t118 = sin(qJ(2));
t144 = qJD(1) * t118;
t92 = -t120 * qJD(2) + t117 * t144;
t107 = t118 * qJDD(1);
t121 = cos(qJ(2));
t142 = qJD(1) * qJD(2);
t139 = t121 * t142;
t97 = t107 + t139;
t67 = -t92 * qJD(3) + t117 * qJDD(2) + t120 * t97;
t104 = t121 * qJD(1) - qJD(3);
t81 = t92 * t104;
t53 = t67 - t81;
t184 = qJ(4) * t53;
t94 = t117 * qJD(2) + t120 * t144;
t73 = t94 * t92;
t106 = t118 * t142;
t141 = t121 * qJDD(1);
t98 = -t106 + t141;
t91 = -qJDD(3) + t98;
t177 = -t73 - t91;
t183 = pkin(3) * t177;
t123 = qJD(2) ^ 2;
t124 = qJD(1) ^ 2;
t167 = t121 * pkin(2);
t133 = -t118 * pkin(6) - t167;
t119 = sin(qJ(1));
t122 = cos(qJ(1));
t132 = t122 * g(1) + t119 * g(2);
t149 = qJDD(1) * pkin(5);
t86 = -t124 * pkin(1) - t132 + t149;
t136 = t124 * t133 + t86;
t166 = t121 * g(3);
t56 = -qJDD(2) * pkin(2) - t123 * pkin(6) + t136 * t118 + t166;
t135 = -t120 * qJDD(2) + t117 * t97;
t66 = -t94 * qJD(3) - t135;
t74 = -t104 * pkin(3) - t94 * qJ(4);
t89 = t92 ^ 2;
t16 = -t66 * pkin(3) - t89 * qJ(4) + t94 * t74 + qJDD(4) + t56;
t182 = t117 * t177;
t181 = t120 * t177;
t130 = t97 + t139;
t131 = -t98 + t106;
t138 = t119 * g(1) - t122 * g(2);
t85 = qJDD(1) * pkin(1) + t124 * pkin(5) + t138;
t47 = t131 * pkin(2) - t130 * pkin(6) - t85;
t168 = t118 * g(3);
t57 = -t123 * pkin(2) + qJDD(2) * pkin(6) + t136 * t121 - t168;
t26 = t117 * t47 + t120 * t57;
t127 = t66 * qJ(4) - 0.2e1 * qJD(4) * t92 + t104 * t74 + t26;
t102 = t104 ^ 2;
t90 = t94 ^ 2;
t69 = -t90 - t102;
t180 = -t127 + (t69 + t89) * pkin(3);
t178 = t67 + t81;
t49 = (qJD(3) + t104) * t94 + t135;
t25 = t117 * t57 - t120 * t47;
t126 = t183 - t25 - t184;
t150 = qJD(4) * t94;
t83 = -0.2e1 * t150;
t9 = t126 + t83;
t175 = pkin(3) * t9;
t68 = -t102 - t89;
t32 = t117 * t68 + t181;
t174 = pkin(2) * t32;
t60 = -t73 + t91;
t156 = t117 * t60;
t36 = t120 * t69 + t156;
t173 = pkin(2) * t36;
t172 = pkin(3) * t53;
t22 = -t117 * t49 - t120 * t53;
t171 = pkin(6) * t22;
t170 = pkin(6) * t32;
t169 = pkin(6) * t36;
t23 = t117 * t53 - t120 * t49;
t59 = -t89 - t90;
t164 = pkin(5) * (t118 * t59 + t121 * t23) - pkin(1) * t22;
t33 = t120 * t68 - t182;
t48 = (qJD(3) - t104) * t94 + t135;
t163 = pkin(5) * (t118 * t48 + t121 * t33) - pkin(1) * t32;
t154 = t120 * t60;
t37 = -t117 * t69 + t154;
t162 = pkin(5) * (t118 * t178 + t121 * t37) - pkin(1) * t36;
t161 = -pkin(2) * t59 + pkin(6) * t23;
t160 = -pkin(2) * t48 + pkin(6) * t33;
t159 = -pkin(2) * t178 + pkin(6) * t37;
t157 = t117 * t56;
t155 = t120 * t56;
t153 = qJ(4) * t117;
t152 = qJ(4) * t120;
t148 = t104 * t117;
t147 = t104 * t120;
t103 = t121 * t124 * t118;
t146 = t118 * (qJDD(2) + t103);
t145 = t121 * (qJDD(2) - t103);
t140 = t121 * t73;
t8 = t117 * t25 + t120 * t26;
t75 = t118 * t86 + t166;
t76 = t121 * t86 - t168;
t137 = t118 * t75 + t121 * t76;
t129 = t117 * t26 - t120 * t25;
t125 = t126 + t183;
t114 = t121 ^ 2;
t113 = t118 ^ 2;
t111 = t114 * t124;
t109 = t113 * t124;
t99 = -0.2e1 * t106 + t141;
t96 = t107 + 0.2e1 * t139;
t84 = 0.2e1 * t150;
t78 = -t90 + t102;
t77 = t89 - t102;
t70 = t90 - t89;
t55 = (t117 * t92 + t120 * t94) * t104;
t42 = t117 * t67 - t94 * t147;
t41 = t120 * t66 - t92 * t148;
t40 = t121 * t91 + t118 * (-t117 * t94 + t120 * t92) * t104;
t39 = t117 * t77 - t154;
t38 = t120 * t78 + t182;
t29 = t118 * (t120 * t67 + t94 * t148) - t140;
t28 = t118 * (-t117 * t66 - t92 * t147) + t140;
t27 = -pkin(3) * t178 + qJ(4) * t60;
t21 = -t117 * t48 + t120 * t178;
t18 = t118 * (t120 * t77 + t156) + t121 * t49;
t17 = t118 * (-t117 * t78 + t181) - t121 * t53;
t13 = t118 * (-t117 * t178 - t120 * t48) - t121 * t70;
t12 = -qJ(4) * t69 + t16;
t10 = -t89 * pkin(3) + t127;
t6 = -pkin(3) * t48 + qJ(4) * t68 - t16;
t5 = -t126 + t84 + t184;
t4 = -qJ(4) * t49 + (-t59 - t89) * pkin(3) + t127;
t3 = -pkin(3) * t16 + qJ(4) * t10;
t2 = t120 * t10 - t117 * t9;
t1 = t117 * t10 + t120 * t9;
t7 = [0, 0, 0, 0, 0, qJDD(1), t138, t132, 0, 0, t130 * t118, t118 * t99 + t121 * t96, t146 + t121 * (-t109 + t123), -t131 * t121, t118 * (t111 - t123) + t145, 0, t121 * t85 + pkin(1) * t99 + pkin(5) * (t121 * (-t111 - t123) - t146), -t118 * t85 - pkin(1) * t96 + pkin(5) * (-t145 - t118 * (-t109 - t123)), pkin(1) * (t109 + t111) + (t113 + t114) * t149 + t137, pkin(1) * t85 + pkin(5) * t137, t29, t13, t17, t28, t18, t40, t118 * (t157 - t170) + t121 * (t25 - t174) + t163, t118 * (t155 - t169) + t121 * (t26 - t173) + t162, t118 * (-t129 - t171) - t22 * t167 + t164, pkin(5) * (t118 * t56 + t121 * t8) + (-pkin(1) + t133) * t129, t29, t13, t17, t28, t18, t40, t118 * (-t117 * t6 - t152 * t177 - t170) + t121 * (-t125 + t84 - t174) + t163, t118 * (-t117 * t27 + t120 * t12 - t169) + t121 * (-t173 - t180) + t162, t118 * (-t117 * t4 + t120 * t5 - t171) + t121 * (-pkin(2) * t22 + t172) + t164, t118 * (-pkin(6) * t1 - t117 * t3 - t9 * t152) + t121 * (-pkin(2) * t1 - t175) - pkin(1) * t1 + pkin(5) * (t118 * t16 + t121 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, t109 - t111, t107, t103, t141, qJDD(2), -t75, -t76, 0, 0, t42, t21, t38, t41, t39, t55, -t155 + t160, t157 + t159, t8 + t161, -pkin(2) * t56 + pkin(6) * t8, t42, t21, t38, t41, t39, t55, t120 * t6 - t153 * t177 + t160, t117 * t12 + t120 * t27 + t159, t117 * t5 + t120 * t4 + t161, -pkin(2) * t16 + pkin(6) * t2 + t120 * t3 - t9 * t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t70, t53, -t73, -t49, -t91, -t25, -t26, 0, 0, t73, t70, t53, -t73, -t49, -t91, t125 + t83, t180, -t172, t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t178, t59, t16;];
tauJ_reg = t7;
