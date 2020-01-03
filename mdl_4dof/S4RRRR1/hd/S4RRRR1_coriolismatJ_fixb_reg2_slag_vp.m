% Calculate inertial parameters regressor of coriolis matrix for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRR1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:17
% EndTime: 2019-12-31 17:22:19
% DurationCPUTime: 1.49s
% Computational Cost: add. (1062->183), mult. (2600->250), div. (0->0), fcn. (1866->6), ass. (0->146)
t105 = sin(qJ(3));
t176 = cos(qJ(2));
t135 = t176 * t105;
t106 = sin(qJ(2));
t108 = cos(qJ(3));
t163 = t106 * t108;
t75 = (t135 + t163) * pkin(1);
t192 = -t75 / 0.2e1;
t104 = sin(qJ(4));
t102 = t104 ^ 2;
t107 = cos(qJ(4));
t103 = t107 ^ 2;
t133 = t103 / 0.2e1 + t102 / 0.2e1;
t134 = t176 * t108;
t164 = t105 * t106;
t77 = (t134 - t164) * pkin(1);
t188 = t133 * t77;
t191 = pkin(3) * t192 + t188 * pkin(7);
t146 = -qJD(2) - qJD(3);
t130 = qJD(1) - t146;
t93 = t103 - t102;
t189 = t130 * t93;
t153 = t102 + t103;
t187 = pkin(3) / 0.2e1;
t145 = t176 * pkin(1);
t127 = t145 + pkin(2);
t83 = t108 * t127;
t70 = pkin(1) * t164 - t83;
t67 = -pkin(3) + t70;
t186 = -t67 / 0.2e1;
t117 = t105 * t127;
t71 = pkin(1) * t163 + t117;
t185 = t71 / 0.2e1;
t184 = -t77 / 0.2e1;
t173 = t108 * pkin(2);
t96 = -pkin(3) - t173;
t183 = -t96 / 0.2e1;
t179 = -t104 / 0.2e1;
t178 = t104 / 0.2e1;
t177 = t107 / 0.2e1;
t175 = pkin(1) * t106;
t174 = pkin(3) * t107;
t39 = t105 * pkin(2);
t66 = t71 * qJD(3);
t72 = t75 * qJD(2);
t172 = -t72 - t66;
t171 = pkin(2) * qJD(2);
t170 = pkin(2) * qJD(3);
t24 = t153 * t70;
t68 = pkin(7) + t71;
t3 = -t68 * t24 + t67 * t71;
t169 = t3 * qJD(1);
t26 = t153 * t77;
t4 = t68 * t26 + t67 * t75;
t168 = t4 * qJD(1);
t167 = t75 * t107;
t8 = t70 * t75 + t71 * t77;
t166 = t8 * qJD(1);
t165 = t96 * t107;
t125 = t39 / 0.2e1 + t185;
t116 = t192 + t125;
t19 = t116 * t107;
t162 = t19 * qJD(1);
t161 = t24 * qJD(1);
t160 = t26 * qJD(1);
t159 = t39 * qJD(1);
t41 = t83 / 0.2e1 + (-t145 / 0.2e1 + pkin(2) / 0.2e1) * t108;
t158 = t41 * qJD(1);
t157 = t70 * qJD(1);
t156 = t71 * qJD(1);
t155 = t75 * qJD(1);
t154 = t77 * qJD(1);
t152 = qJD(1) * t104;
t151 = qJD(1) * t107;
t150 = qJD(2) * t104;
t149 = qJD(3) * t104;
t148 = t104 * qJD(4);
t101 = t107 * qJD(4);
t147 = -qJD(1) - qJD(2);
t144 = t105 * t170;
t143 = t105 * t171;
t142 = -t173 / 0.2e1;
t141 = t67 * t152;
t140 = t67 * t151;
t139 = t71 * t152;
t138 = t75 * t152;
t132 = t176 * qJD(1);
t131 = t176 * qJD(2);
t129 = pkin(2) * t146;
t128 = t153 * t108;
t126 = t70 / 0.2e1 + t187 + t186;
t124 = t184 + t183 + t186;
t123 = t133 * t70;
t122 = t105 * t129;
t121 = t133 * t108;
t95 = pkin(7) + t39;
t109 = (t67 * t105 / 0.2e1 + t68 * t121) * pkin(2) - t95 * t123 + t96 * t185;
t2 = t109 - t191;
t27 = (t105 * t96 + t95 * t128) * pkin(2);
t120 = -t2 * qJD(1) - t27 * qJD(2);
t6 = t153 * (t173 / 0.2e1 - t70 / 0.2e1 + t184);
t76 = pkin(2) * t128;
t119 = -t6 * qJD(1) - t76 * qJD(2);
t118 = t142 + t187 + t183;
t9 = t124 * t104;
t114 = t9 * qJD(1) - t96 * t150;
t10 = t124 * t107;
t113 = t10 * qJD(1) - qJD(2) * t165;
t18 = t116 * t104;
t112 = -t18 * qJD(1) - t104 * t143;
t13 = t126 * t104;
t38 = t118 * t104;
t111 = pkin(3) * t149 + t13 * qJD(1) + t38 * qJD(2);
t14 = t126 * t107;
t40 = t118 * t107;
t110 = t14 * qJD(1) + t40 * qJD(2) + qJD(3) * t174;
t100 = -t174 / 0.2e1;
t99 = pkin(3) * t179;
t94 = t104 * t101;
t92 = t104 * t144;
t82 = t93 * qJD(4);
t81 = t165 / 0.2e1;
t80 = t96 * t178;
t74 = t77 * qJD(2);
t73 = t76 * qJD(3);
t65 = t70 * qJD(3);
t64 = t75 * t150;
t55 = t71 * t149;
t48 = t67 * t177;
t47 = t67 * t178;
t44 = t130 * t107 * t104;
t43 = t107 * t142 + t100 + t81;
t42 = t104 * t142 + t80 + t99;
t29 = t142 - t83 / 0.2e1 + (t164 - t134 / 0.2e1) * pkin(1);
t28 = -t39 / 0.2e1 - t117 / 0.2e1 + (-t163 - t135 / 0.2e1) * pkin(1);
t25 = t26 * qJD(2);
t23 = t24 * qJD(3);
t20 = -t167 / 0.2e1 - t125 * t107;
t17 = t125 * t104 + t75 * t178;
t16 = t70 * t177 + t100 + t48;
t15 = t70 * t178 + t47 + t99;
t12 = t107 * t184 + t48 + t81;
t11 = t77 * t179 + t47 + t80;
t5 = pkin(2) * t121 - t123 + t188;
t1 = t109 + t191;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t175, -pkin(1) * t131, 0, 0, 0, 0, 0, 0, 0, 0, t172, -t74 + t65, 0, t8 * qJD(2), t94, t82, 0, -t94, 0, 0, t172 * t107 + t67 * t148, t101 * t67 + t55 + t64, t25 - t23, t4 * qJD(2) + t3 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147 * t175, (-t132 - t131) * pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, t28 * qJD(3) - t155 - t72, t29 * qJD(3) - t154 - t74, 0, t166 + (t105 * t77 - t108 * t75) * t171, t94, t82, 0, -t94, 0, 0, t20 * qJD(3) + t11 * qJD(4) + t147 * t167, t17 * qJD(3) + t12 * qJD(4) + t138 + t64, t5 * qJD(3) + t160 + t25, t168 + (t26 * t95 + t75 * t96) * qJD(2) + t1 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * qJD(2) - t156 - t66, t29 * qJD(2) + t157 + t65, 0, 0, t94, t82, 0, -t94, 0, 0, t20 * qJD(2) + t15 * qJD(4) + (-qJD(1) - qJD(3)) * t71 * t107, t17 * qJD(2) + t16 * qJD(4) + t139 + t55, t5 * qJD(2) - t161 - t23, t169 + t1 * qJD(2) + (-t71 * pkin(3) - pkin(7) * t24) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t189, t101, -t44, -t148, 0, t11 * qJD(2) + t15 * qJD(3) - t101 * t68 + t141, t12 * qJD(2) + t16 * qJD(3) + t148 * t68 + t140, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1) * t175, pkin(1) * t132, 0, 0, 0, 0, 0, 0, 0, 0, -t39 * qJD(3) + t155, -t41 * qJD(3) + t154, 0, -t166, t94, t82, 0, -t94, 0, 0, -t19 * qJD(3) - t9 * qJD(4) + t151 * t75, t18 * qJD(3) - t10 * qJD(4) - t138, t6 * qJD(3) - t160, t2 * qJD(3) - t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, -t108 * t170, 0, 0, t94, t82, 0, -t94, 0, 0, -t107 * t144 + t148 * t96, t101 * t96 + t92, t73, t27 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122 - t159, t108 * t129 - t158, 0, 0, t94, t82, 0, -t94, 0, 0, t42 * qJD(4) + t107 * t122 - t162, t43 * qJD(4) - t112 + t92, -t119 + t73, (-pkin(3) * t105 + pkin(7) * t128) * t170 - t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t189, t101, -t44, -t148, 0, t42 * qJD(3) - t101 * t95 - t114, t43 * qJD(3) + t148 * t95 - t113, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * qJD(2) + t156, t41 * qJD(2) - t157, 0, 0, t94, t82, 0, -t94, 0, 0, t19 * qJD(2) - t13 * qJD(4) + t151 * t71, -t18 * qJD(2) - t14 * qJD(4) - t139, -t6 * qJD(2) + t161, -t2 * qJD(2) - t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143 + t159, t108 * t171 + t158, 0, 0, t94, t82, 0, -t94, 0, 0, -t38 * qJD(4) + t107 * t143 + t162, -t40 * qJD(4) + t112, t119, t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t82, 0, -t94, 0, 0, -pkin(3) * t148, -pkin(3) * t101, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t189, t101, -t44, -t148, 0, -pkin(7) * t101 - t111, pkin(7) * t148 - t110, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t189, 0, t44, 0, 0, t9 * qJD(2) + t13 * qJD(3) - t141, t10 * qJD(2) + t14 * qJD(3) - t140, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t189, 0, t44, 0, 0, t38 * qJD(3) + t114, t40 * qJD(3) + t113, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t189, 0, t44, 0, 0, t111, t110, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
