% Calculate inertial parameters regressor of coriolis matrix for
% S4RRRP4
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
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:44
% EndTime: 2019-12-31 17:15:47
% DurationCPUTime: 1.16s
% Computational Cost: add. (1430->135), mult. (2910->169), div. (0->0), fcn. (2928->4), ass. (0->115)
t127 = qJD(2) + qJD(3);
t163 = cos(qJ(3));
t125 = t163 * pkin(2);
t106 = sin(qJ(3));
t107 = sin(qJ(2));
t108 = cos(qJ(2));
t92 = t106 * t107 - t163 * t108;
t149 = t92 * qJ(4);
t168 = -pkin(6) - pkin(5);
t98 = t168 * t107;
t157 = t106 * t98;
t99 = t168 * t108;
t97 = t163 * t99;
t68 = -t97 + t157;
t40 = t68 - t149;
t172 = -t106 * t99 - t163 * t98;
t94 = t106 * t108 + t163 * t107;
t175 = t94 * qJ(4) + t172;
t177 = t127 * t175;
t173 = t127 * t94;
t176 = t92 * t173;
t174 = t127 * t172;
t171 = t94 ^ 2;
t170 = -pkin(3) / 0.2e1;
t169 = t97 / 0.2e1;
t167 = t40 * pkin(3);
t166 = t92 * pkin(3);
t165 = t94 * pkin(3);
t103 = t125 + pkin(3);
t164 = -t103 / 0.2e1;
t162 = pkin(2) * t106;
t161 = t107 * pkin(2);
t104 = -t108 * pkin(2) - pkin(1);
t71 = t104 + t166;
t51 = t71 * t94;
t159 = pkin(3) * qJD(3);
t158 = qJD(2) * pkin(2);
t72 = t161 + t165;
t5 = t71 * t72;
t153 = t5 * qJD(1);
t6 = pkin(3) * t51;
t152 = t6 * qJD(1);
t9 = t104 * t161;
t150 = t9 * qJD(1);
t10 = t175 * t94 - t40 * t92;
t148 = t10 * qJD(1);
t13 = t72 * t92 + t51;
t147 = t13 * qJD(1);
t50 = t71 * t92;
t14 = t72 * t94 - t50;
t146 = t14 * qJD(1);
t17 = -t92 * t165 - t51;
t145 = t17 * qJD(1);
t18 = -t171 * pkin(3) + t50;
t144 = t18 * qJD(1);
t118 = -t106 * t92 / 0.2e1;
t20 = (t164 + t170) * t94 + (t118 - t107 / 0.2e1) * pkin(2);
t143 = t20 * qJD(1);
t110 = -t125 / 0.2e1 + t103 / 0.2e1;
t26 = (t170 + t110) * t92;
t142 = t26 * qJD(1);
t141 = t40 * qJD(3);
t91 = t92 ^ 2;
t41 = t91 - t171;
t140 = t41 * qJD(1);
t45 = t104 * t94 + t92 * t161;
t139 = t45 * qJD(1);
t46 = -t104 * t92 + t94 * t161;
t138 = t46 * qJD(1);
t60 = t91 + t171;
t137 = t60 * qJD(1);
t64 = t169 - t97 / 0.2e1;
t136 = t64 * qJD(1);
t135 = t92 * qJD(1);
t134 = t94 * qJD(1);
t84 = t94 * qJD(4);
t133 = qJD(1) * t104;
t132 = qJD(1) * t108;
t131 = qJD(3) * t104;
t100 = -t107 ^ 2 + t108 ^ 2;
t130 = t100 * qJD(1);
t129 = t107 * qJD(2);
t128 = t108 * qJD(2);
t126 = pkin(3) * t134;
t124 = pkin(1) * t107 * qJD(1);
t123 = pkin(1) * t132;
t122 = qJD(3) * t162;
t121 = t92 * t134;
t120 = t92 * t133;
t119 = t94 * t133;
t117 = t107 * t128;
t116 = t163 * qJD(2);
t115 = t163 * qJD(3);
t113 = pkin(2) * t115;
t109 = (t125 / 0.2e1 + t164) * t40;
t2 = t167 / 0.2e1 + t109;
t76 = (t125 - t103) * t162;
t111 = -t2 * qJD(1) - t76 * qJD(2);
t43 = 0.2e1 * t169 - t157;
t101 = t107 * t132;
t80 = t92 * qJD(4);
t77 = t94 * t162;
t59 = t64 * qJD(2);
t58 = t64 * qJD(3);
t54 = t127 * t92;
t49 = t106 * t158 - t136;
t48 = pkin(2) * t116;
t35 = -t127 * t162 + t136;
t34 = (-t116 - t115) * pkin(2);
t25 = t166 / 0.2e1 + t110 * t92;
t24 = t43 + t149;
t19 = pkin(2) * t118 + t94 * t164 + t161 / 0.2e1 + t165 / 0.2e1;
t12 = t127 * t41;
t1 = -t167 / 0.2e1 + t109;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, t100 * qJD(2), 0, -t117, 0, 0, -pkin(1) * t129, -pkin(1) * t128, 0, 0, -t176, t12, 0, t176, 0, 0, t45 * qJD(2) + t94 * t131, t46 * qJD(2) - t92 * t131, 0, t9 * qJD(2), -t176, t12, 0, t176, 0, 0, t13 * qJD(2) - t17 * qJD(3), t14 * qJD(2) - t18 * qJD(3), t60 * qJD(4), t5 * qJD(2) + t6 * qJD(3) + t10 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t130, t128, -t101, -t129, 0, -pkin(5) * t128 - t124, pkin(5) * t129 - t123, 0, 0, -t121, t140, -t54, t121, -t173, 0, -qJD(2) * t68 + t43 * qJD(3) + t139, t138 + t174, (t92 * t125 - t77) * qJD(2), t150 + (-t106 * t172 - t163 * t68) * t158, -t121, t140, -t54, t121, -t173, 0, -qJD(2) * t40 + t24 * qJD(3) + t147, t146 + t177, (t103 * t92 - t77) * qJD(2) + t25 * qJD(3), t153 + (-t40 * t103 - t162 * t175) * qJD(2) + t1 * qJD(3) + t19 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t140, -t54, t121, -t173, 0, t43 * qJD(2) - t68 * qJD(3) + t119, -t120 + t174, 0, 0, -t121, t140, -t54, t121, -t173, 0, t24 * qJD(2) - t141 - t145, -t144 + t177, t25 * qJD(2) + t92 * t159, -pkin(3) * t141 + t1 * qJD(2) + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, t19 * qJD(2) + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t130, 0, t101, 0, 0, t124, t123, 0, 0, t121, -t140, 0, -t121, 0, 0, t58 - t139, -t138, 0, -t150, t121, -t140, 0, -t121, 0, 0, t58 - t84 - t147, t80 - t146, t26 * qJD(3), t2 * qJD(3) + t20 * qJD(4) - t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t113, 0, 0, 0, 0, 0, 0, 0, 0, -t122, -t113, 0, t76 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t34, 0, 0, 0, 0, 0, 0, 0, 0, t35, t34, t142, -pkin(3) * t122 - t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, t135, 0, t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, -t140, 0, -t121, 0, 0, -t59 - t119, t120, 0, 0, t121, -t140, 0, -t121, 0, 0, -t59 - t84 + t145, t80 + t144, -t26 * qJD(2), -pkin(3) * t84 - t2 * qJD(2) - t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t48, 0, 0, 0, 0, 0, 0, 0, 0, t49, t48, -t142, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, t135, 0, -t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, -t54, -t137, -t20 * qJD(2) + t94 * t159 - t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t135, 0, -t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t135, 0, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
