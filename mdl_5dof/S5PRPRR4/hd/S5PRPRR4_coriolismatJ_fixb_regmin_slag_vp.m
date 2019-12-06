% Calculate minimal parameter regressor of coriolis matrix for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRPRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:38
% EndTime: 2019-12-05 15:51:41
% DurationCPUTime: 0.98s
% Computational Cost: add. (757->131), mult. (2270->243), div. (0->0), fcn. (2411->10), ass. (0->120)
t93 = sin(qJ(4));
t138 = t93 * qJD(2);
t95 = cos(qJ(5));
t131 = t95 * t138;
t92 = sin(qJ(5));
t86 = t92 ^ 2;
t88 = t95 ^ 2;
t77 = t88 - t86;
t167 = qJD(4) * t77 - 0.2e1 * t92 * t131;
t96 = cos(qJ(4));
t162 = t96 * pkin(8);
t163 = t93 * pkin(4);
t72 = -t162 + t163;
t166 = -t72 / 0.2e1;
t165 = -t93 / 0.2e1;
t164 = -t96 / 0.2e1;
t161 = cos(qJ(2));
t150 = sin(pkin(5));
t151 = cos(pkin(10));
t110 = t151 * t150;
t94 = sin(qJ(2));
t118 = t94 * t150;
t90 = sin(pkin(10));
t57 = -t110 * t161 + t118 * t90;
t23 = t57 * t96;
t113 = t150 * t161;
t58 = t110 * t94 + t113 * t90;
t160 = t58 * t92;
t159 = t58 * t95;
t87 = t93 ^ 2;
t158 = t92 * t87;
t157 = t92 * t93;
t156 = t92 * t96;
t155 = t93 * t95;
t154 = t95 * t72;
t153 = t95 * t87;
t152 = t95 * t96;
t89 = t96 ^ 2;
t78 = t89 - t87;
t149 = qJD(4) * t92;
t148 = qJD(4) * t95;
t147 = qJD(5) * t92;
t146 = qJD(5) * t95;
t145 = qJD(5) * t96;
t82 = pkin(2) * t90 + pkin(7);
t133 = t82 * t156;
t114 = -t96 * pkin(4) - t93 * pkin(8);
t83 = -pkin(2) * t151 - pkin(3);
t61 = t114 + t83;
t41 = -t61 * t95 + t133;
t70 = t82 * t157;
t10 = t41 * t93 + (-t70 + t154) * t96;
t144 = t10 * qJD(2);
t132 = t82 * t152;
t42 = t61 * t92 + t132;
t11 = t72 * t156 + (-t42 + t132) * t93;
t143 = t11 * qJD(2);
t91 = cos(pkin(5));
t44 = t58 * t96 + t91 * t93;
t142 = t44 * qJD(4);
t67 = t78 * t92;
t141 = t67 * qJD(2);
t68 = t89 * t95 - t153;
t140 = t68 * qJD(2);
t139 = t78 * qJD(2);
t137 = t93 * qJD(4);
t136 = t93 * qJD(5);
t135 = t96 * qJD(2);
t134 = t96 * qJD(4);
t130 = t92 * t145;
t129 = t95 * t145;
t128 = t83 * t138;
t127 = t83 * t135;
t126 = t92 * t146;
t125 = t92 * t148;
t124 = t93 * t134;
t123 = t93 * t135;
t122 = t95 * t137;
t121 = t23 / 0.2e1;
t120 = -t157 / 0.2e1;
t119 = t155 / 0.2e1;
t116 = t92 * t122;
t115 = -qJD(5) + t135;
t112 = t44 * t95 + t57 * t92;
t111 = t44 * t92 - t57 * t95;
t24 = -t158 * t82 - t41 * t96;
t43 = t58 * t93 - t91 * t96;
t97 = t58 / 0.2e1 + t44 * t164 + t43 * t165;
t6 = t97 * t92;
t109 = -qJD(1) * t6 - qJD(2) * t24;
t25 = -t153 * t82 - t42 * t96;
t5 = t97 * t95;
t108 = qJD(1) * t5 + qJD(2) * t25;
t107 = t115 * t93;
t106 = t162 / 0.2e1 - t163 / 0.2e1;
t101 = t166 + t106;
t45 = t101 * t92;
t105 = pkin(4) * t148 + qJD(2) * t45;
t46 = t101 * t95;
t104 = pkin(4) * t149 - qJD(2) * t46;
t103 = t95 * t107;
t62 = (t86 / 0.2e1 - t88 / 0.2e1) * t93;
t102 = -qJD(2) * t62 + t125;
t100 = qJD(2) * t153 * t92 + qJD(4) * t62;
t66 = t77 * t87;
t99 = qJD(2) * t66 + 0.2e1 * t116;
t84 = t137 / 0.2e1;
t81 = t92 * t137;
t65 = (t135 - qJD(5) / 0.2e1) * t93;
t59 = t62 * qJD(5);
t27 = t70 + t154 / 0.2e1 + t106 * t95;
t26 = t82 * t155 + (-t106 + t166) * t92;
t21 = t57 * t93;
t15 = t43 * t95;
t13 = t43 * t92;
t8 = t112 * t96 / 0.2e1 + t43 * t119 + t92 * t121 + t159 / 0.2e1;
t7 = t111 * t164 + t43 * t120 + t95 * t121 - t160 / 0.2e1;
t4 = t112 * t165 + t119 * t44 + t120 * t57;
t3 = t111 * t165 + t44 * t157 / 0.2e1 + t57 * t119;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -qJD(2) * t118, -qJD(2) * t113, (-t151 * t58 - t57 * t90) * qJD(2) * pkin(2), 0, 0, 0, 0, 0, qJD(4) * t21 - t135 * t58, qJD(4) * t23 + t138 * t58, 0, 0, 0, 0, 0, (-(t156 * t57 + t159) * t96 - t57 * t158) * qJD(2) + t3 * qJD(4) + t8 * qJD(5), ((-t152 * t57 + t160) * t96 - t57 * t153) * qJD(2) + t4 * qJD(4) + t7 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t21 - t142, qJD(2) * t23 + qJD(4) * t43, 0, 0, 0, 0, 0, qJD(2) * t3 + qJD(5) * t13 - t142 * t95, qJD(2) * t4 + qJD(5) * t15 + t142 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(2) + t13 * qJD(4) - qJD(5) * t112, t7 * qJD(2) + t15 * qJD(4) + qJD(5) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * qJD(5), t6 * qJD(5); 0, 0, 0, 0, 0, t124, t78 * qJD(4), 0, 0, 0, t83 * t137, t83 * t134, t124 * t88 - t126 * t87, -qJD(5) * t66 - 0.2e1 * t116 * t96, -qJD(4) * t68 + t130 * t93, qJD(4) * t67 + t129 * t93, -t124, -qJD(4) * t10 - qJD(5) * t25, qJD(4) * t11 + qJD(5) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t123, t139, t134, -t137, 0, -t134 * t82 + t128, t137 * t82 + t127, -t59 + (t138 * t88 + t125) * t96, -0.2e1 * t93 * t126 + t167 * t96, t81 - t140, t122 + t141, -t65, -t144 + (t114 * t92 - t132) * qJD(4) + t27 * qJD(5), t143 + (t114 * t95 + t133) * qJD(4) + t26 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t99, t92 * t107, t103, t84, qJD(4) * t27 - qJD(5) * t42 - t108, qJD(4) * t26 + qJD(5) * t41 - t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, -t134, 0, 0, 0, 0, 0, -t122 - t130, t81 - t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134 * t92 - t136 * t95, -t134 * t95 + t136 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t123, -t139, 0, 0, 0, -t128, -t127, -t123 * t88 - t59, 0.2e1 * t92 * t103, -t129 + t140, t130 - t141, t65, qJD(5) * t46 + t144, -qJD(5) * t45 - t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, t77 * qJD(5), 0, 0, 0, -pkin(4) * t147, -pkin(4) * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t167, -t115 * t95, t115 * t92, -t138 / 0.2e1, -pkin(8) * t146 - t104, pkin(8) * t147 - t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(2), -t6 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t99, (-t138 * t92 + t148) * t96, (-t131 - t149) * t96, t84, -qJD(4) * t46 + t108, qJD(4) * t45 + t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t167, t95 * t135, -t92 * t135, t138 / 0.2e1, t104, t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
