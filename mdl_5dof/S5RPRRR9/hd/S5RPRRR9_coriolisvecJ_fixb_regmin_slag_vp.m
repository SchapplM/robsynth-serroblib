% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR9
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
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:59
% EndTime: 2019-12-31 19:08:06
% DurationCPUTime: 1.87s
% Computational Cost: add. (2933->227), mult. (7928->318), div. (0->0), fcn. (6378->8), ass. (0->135)
t103 = cos(qJ(5));
t138 = qJD(5) * t103;
t101 = sin(qJ(4));
t104 = cos(qJ(4));
t105 = cos(qJ(3));
t99 = cos(pkin(9));
t146 = t105 * t99;
t102 = sin(qJ(3));
t98 = sin(pkin(9));
t150 = t102 * t98;
t116 = -t146 + t150;
t76 = t116 * qJD(1);
t82 = t102 * t99 + t105 * t98;
t77 = t82 * qJD(1);
t57 = t101 * t77 + t104 * t76;
t184 = t103 * t57;
t192 = t138 + t184;
t97 = qJD(3) + qJD(4);
t160 = t57 * t97;
t140 = qJD(4) * t104;
t141 = qJD(4) * t101;
t135 = qJD(1) * t150;
t142 = qJD(3) * t105;
t88 = t99 * qJD(1) * t142;
t72 = -qJD(3) * t135 + t88;
t79 = t82 * qJD(3);
t73 = qJD(1) * t79;
t32 = -t101 * t73 + t104 * t72 - t76 * t140 - t77 * t141;
t191 = t32 + t160;
t181 = -qJD(5) - t57;
t190 = qJD(5) + t181;
t100 = sin(qJ(5));
t119 = -t101 * t76 + t104 * t77;
t139 = qJD(5) * t100;
t14 = t103 * t32 - t119 * t139 + t97 * t138;
t51 = t100 * t97 + t103 * t119;
t15 = qJD(5) * t51 + t100 * t32;
t49 = t100 * t119 - t103 * t97;
t189 = -t100 * t15 + t14 * t103 - t192 * t49;
t12 = t14 * t100;
t188 = t192 * t51 + t12;
t33 = qJD(4) * t119 + t101 * t72 + t104 * t73;
t29 = t100 * t33;
t52 = t181 * t138;
t157 = t29 - t52;
t162 = t51 * t119;
t187 = -t181 * t184 + t157 - t162;
t159 = pkin(6) + qJ(2);
t86 = t159 * t98;
t83 = qJD(1) * t86;
t87 = t159 * t99;
t84 = qJD(1) * t87;
t118 = t102 * t83 - t105 * t84;
t48 = -pkin(7) * t76 - t118;
t152 = t101 * t48;
t173 = -t102 * t84 - t105 * t83;
t47 = -pkin(7) * t77 + t173;
t46 = qJD(3) * pkin(3) + t47;
t20 = t104 * t46 - t152;
t18 = -pkin(4) * t97 - t20;
t186 = t18 * t57;
t92 = -pkin(2) * t99 - pkin(1);
t85 = t92 * qJD(1) + qJD(2);
t63 = pkin(3) * t76 + t85;
t185 = t57 * t63;
t183 = t119 * t57;
t182 = t100 * t181;
t161 = t119 * t97;
t180 = -t33 + t161;
t178 = t119 ^ 2 - t57 ^ 2;
t35 = pkin(4) * t119 + pkin(8) * t57;
t163 = t49 * t119;
t175 = t181 * t119;
t31 = t103 * t33;
t174 = -t139 * t181 - t31;
t147 = t104 * t48;
t21 = t101 * t46 + t147;
t19 = pkin(8) * t97 + t21;
t22 = pkin(4) * t57 - pkin(8) * t119 + t63;
t120 = t100 * t19 - t103 * t22;
t172 = t119 * t120 + t18 * t139;
t39 = -pkin(7) * t73 - qJD(2) * t76 + t173 * qJD(3);
t113 = t82 * qJD(2);
t111 = qJD(1) * t113;
t40 = -pkin(7) * t72 + qJD(3) * t118 - t111;
t128 = t101 * t39 - t104 * t40;
t3 = qJD(4) * t21 + t128;
t5 = t100 * t22 + t103 * t19;
t171 = t3 * t100 + t5 * t119 + t18 * t138;
t170 = -t119 * t63 - t128;
t127 = t101 * t40 - t48 * t141;
t2 = (qJD(4) * t46 + t39) * t104 + t127;
t53 = -pkin(7) * t82 - t102 * t87 - t105 * t86;
t117 = t102 * t86 - t105 * t87;
t54 = -pkin(7) * t116 - t117;
t27 = t101 * t53 + t104 * t54;
t61 = t101 * t82 + t104 * t116;
t62 = -t101 * t116 + t104 * t82;
t67 = pkin(3) * t116 + t92;
t28 = pkin(4) * t61 - pkin(8) * t62 + t67;
t78 = t116 * qJD(3);
t37 = -qJD(4) * t61 - t101 * t79 - t104 * t78;
t26 = t101 * t54 - t104 * t53;
t110 = -t86 * t142 + qJD(2) * t146 + (-qJD(2) * t98 - qJD(3) * t87) * t102;
t41 = -pkin(7) * t79 + t110;
t107 = qJD(3) * t117 - t113;
t42 = pkin(7) * t78 + t107;
t6 = -qJD(4) * t26 + t101 * t42 + t104 * t41;
t169 = (qJD(5) * t28 + t6) * t181 - (qJD(5) * t22 + t2) * t61 + t18 * t37 - t27 * t33 + t3 * t62;
t168 = pkin(3) * t77;
t166 = t18 * t62;
t165 = t28 * t33;
t164 = t33 * t62;
t156 = t98 ^ 2 + t99 ^ 2;
t154 = t100 * t51;
t137 = qJD(1) * qJD(2);
t133 = t62 * t139;
t132 = -pkin(3) * t97 - t46;
t131 = t156 * qJD(1) ^ 2;
t93 = pkin(3) * t101 + pkin(8);
t124 = qJD(5) * t93 + t168 + t35;
t23 = t101 * t47 + t147;
t122 = pkin(3) * t141 - t23;
t121 = -t181 * t37 + t164;
t115 = 0.2e1 * t156 * t137;
t114 = t182 * t57 - t174;
t24 = t104 * t47 - t152;
t108 = t186 - t93 * t33 - (-pkin(3) * t140 + t24) * t181;
t94 = -pkin(3) * t104 - pkin(4);
t38 = qJD(4) * t62 - t101 * t78 + t104 * t79;
t10 = pkin(3) * t79 + pkin(4) * t38 - pkin(8) * t37;
t9 = pkin(3) * t73 + pkin(4) * t33 - pkin(8) * t32;
t8 = t103 * t9;
t7 = qJD(4) * t27 + t101 * t41 - t104 * t42;
t1 = [0, 0, 0, 0, 0, t115, qJ(2) * t115, t72 * t82 - t77 * t78, -t116 * t72 - t73 * t82 + t76 * t78 - t77 * t79, -t78 * qJD(3), -t79 * qJD(3), 0, qJD(3) * t107 + t92 * t73 + t85 * t79, -qJD(3) * t110 + t92 * t72 - t85 * t78, t119 * t37 + t32 * t62, -t119 * t38 - t32 * t61 - t37 * t57 - t164, t37 * t97, -t38 * t97, 0, t33 * t67 + t38 * t63 - t7 * t97 + (t57 * t79 + t61 * t73) * pkin(3), t32 * t67 + t37 * t63 - t6 * t97 + (t119 * t79 + t62 * t73) * pkin(3), -t51 * t133 + (t14 * t62 + t37 * t51) * t103, (-t103 * t49 - t154) * t37 + (-t12 - t103 * t15 + (t100 * t49 - t103 * t51) * qJD(5)) * t62, t103 * t121 + t133 * t181 + t14 * t61 + t38 * t51, -t100 * t121 - t15 * t61 - t38 * t49 + t52 * t62, -t181 * t38 + t33 * t61, t26 * t15 - t120 * t38 + t7 * t49 + t8 * t61 + (-t10 * t181 + t165 + (t181 * t27 - t19 * t61 + t166) * qJD(5)) * t103 + t169 * t100, t26 * t14 - t5 * t38 + t7 * t51 + ((-qJD(5) * t27 + t10) * t181 - t165 - (-qJD(5) * t19 + t9) * t61 - qJD(5) * t166) * t100 + t169 * t103; 0, 0, 0, 0, 0, -t131, -qJ(2) * t131, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t77, t88 + (-t76 - t135) * qJD(3), 0, 0, 0, 0, 0, t33 + t161, t32 - t160, 0, 0, 0, 0, 0, t114 - t163, -t103 * t181 ^ 2 - t162 - t29; 0, 0, 0, 0, 0, 0, 0, t77 * t76, -t76 ^ 2 + t77 ^ 2, t88 + (t76 - t135) * qJD(3), 0, 0, -t85 * t77 - t111, t116 * t137 + t85 * t76, t183, t178, t191, t180, 0, -t57 * t168 + t23 * t97 + (t101 * t132 - t147) * qJD(4) + t170, -t119 * t168 + t24 * t97 + t185 + (qJD(4) * t132 - t39) * t104 - t127, t188, t154 * t181 + t189, t187, t114 + t163, t175, t94 * t15 + t122 * t49 + (t124 * t181 - t3) * t103 + t108 * t100 + t172, t103 * t108 + t122 * t51 - t124 * t182 + t94 * t14 + t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, t178, t191, t180, 0, (-qJD(4) + t97) * t21 + t170, t20 * t97 + t185 - t2, t188, t182 * t51 + t189, t187, -t181 * t182 + t163 + t31, t175, -pkin(4) * t15 - t3 * t103 + (-t100 * t20 + t103 * t35) * t181 - t21 * t49 + t100 * t186 - t157 * pkin(8) + t172, -pkin(4) * t14 - (t100 * t35 + t103 * t20) * t181 - t21 * t51 + t18 * t184 + t174 * pkin(8) + t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51 * t49, -t49 ^ 2 + t51 ^ 2, -t181 * t49 + t14, -t181 * t51 - t15, t33, -t100 * t2 - t18 * t51 - t190 * t5 + t8, -t100 * t9 - t103 * t2 + t190 * t120 + t18 * t49;];
tauc_reg = t1;
