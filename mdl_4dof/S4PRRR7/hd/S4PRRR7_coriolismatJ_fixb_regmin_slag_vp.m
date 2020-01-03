% Calculate minimal parameter regressor of coriolis matrix for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:49
% EndTime: 2019-12-31 16:36:52
% DurationCPUTime: 0.85s
% Computational Cost: add. (462->120), mult. (1462->230), div. (0->0), fcn. (1403->8), ass. (0->114)
t78 = sin(qJ(3));
t121 = t78 * qJD(2);
t80 = cos(qJ(4));
t114 = t80 * t121;
t77 = sin(qJ(4));
t71 = t77 ^ 2;
t73 = t80 ^ 2;
t63 = t73 - t71;
t153 = t63 * qJD(3) - 0.2e1 * t77 * t114;
t81 = cos(qJ(3));
t148 = t81 * pkin(7);
t149 = t78 * pkin(3);
t58 = -t148 + t149;
t152 = -t58 / 0.2e1;
t151 = -t78 / 0.2e1;
t150 = -t81 / 0.2e1;
t75 = sin(pkin(4));
t79 = sin(qJ(2));
t147 = t75 * t79;
t82 = cos(qJ(2));
t146 = t75 * t82;
t145 = t77 * t78;
t144 = t77 * t79;
t143 = t77 * t81;
t142 = t78 * t80;
t141 = t79 * t80;
t140 = t80 * t58;
t72 = t78 ^ 2;
t139 = t80 * t72;
t138 = t81 * t82;
t137 = t82 * t72;
t74 = t81 ^ 2;
t64 = t74 - t72;
t119 = pkin(6) * t143;
t98 = -t81 * pkin(3) - t78 * pkin(7);
t57 = -pkin(2) + t98;
t34 = -t80 * t57 + t119;
t68 = pkin(6) * t145;
t9 = t34 * t78 + (-t68 + t140) * t81;
t136 = t9 * qJD(2);
t135 = qJD(2) * t75;
t134 = qJD(3) * t77;
t133 = qJD(3) * t78;
t132 = qJD(3) * t80;
t131 = qJD(3) * t81;
t130 = qJD(3) * t82;
t129 = qJD(4) * t77;
t128 = qJD(4) * t80;
t127 = qJD(4) * t81;
t118 = t80 * t81 * pkin(6);
t35 = t77 * t57 + t118;
t10 = t58 * t143 + (-t35 + t118) * t78;
t126 = t10 * qJD(2);
t76 = cos(pkin(4));
t42 = t81 * t147 + t76 * t78;
t125 = t42 * qJD(3);
t52 = t64 * t77;
t124 = t52 * qJD(2);
t53 = t80 * t74 - t139;
t123 = t53 * qJD(2);
t122 = t64 * qJD(2);
t120 = t81 * qJD(2);
t117 = t80 * t146;
t116 = pkin(2) * t121;
t115 = pkin(2) * t120;
t113 = t77 * t127;
t112 = t80 * t127;
t111 = t82 * t135;
t110 = t77 * t128;
t109 = t77 * t132;
t108 = t78 * t131;
t107 = t78 * t120;
t106 = t78 * t132;
t41 = t78 * t147 - t76 * t81;
t105 = t41 * t151;
t104 = t145 / 0.2e1;
t103 = t142 / 0.2e1;
t102 = -t138 / 0.2e1;
t100 = t77 * t106;
t99 = -qJD(4) + t120;
t21 = -t72 * pkin(6) * t77 - t34 * t81;
t83 = t147 / 0.2e1 + t42 * t150 + t105;
t6 = t83 * t77;
t97 = -t6 * qJD(1) - t21 * qJD(2);
t22 = -pkin(6) * t139 - t35 * t81;
t5 = t83 * t80;
t96 = t5 * qJD(1) + t22 * qJD(2);
t95 = t99 * t78;
t94 = t148 / 0.2e1 - t149 / 0.2e1;
t93 = t42 * t77 + t117;
t92 = t77 * t146 - t42 * t80;
t87 = t152 + t94;
t23 = t87 * t77;
t91 = pkin(3) * t132 + t23 * qJD(2);
t24 = t87 * t80;
t90 = pkin(3) * t134 - t24 * qJD(2);
t89 = t80 * t95;
t43 = (t71 / 0.2e1 - t73 / 0.2e1) * t78;
t88 = -t43 * qJD(2) + t109;
t86 = t77 * qJD(2) * t139 + t43 * qJD(3);
t51 = t63 * t72;
t85 = t51 * qJD(2) + 0.2e1 * t100;
t69 = t133 / 0.2e1;
t46 = (t120 - qJD(4) / 0.2e1) * t78;
t40 = t43 * qJD(4);
t16 = t68 + t140 / 0.2e1 + t94 * t80;
t15 = pkin(6) * t142 + (t152 - t94) * t77;
t14 = t41 * t80;
t12 = t41 * t77;
t8 = t92 * t150 + t41 * t103 + (t77 * t102 + t141 / 0.2e1) * t75;
t7 = t93 * t150 + t77 * t105 + (t80 * t102 - t144 / 0.2e1) * t75;
t4 = t92 * t78 / 0.2e1 + t42 * t103 + t104 * t146;
t3 = t42 * t104 + (t117 + t93) * t151;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t79 * t135, -t111, 0, 0, 0, 0, 0, (-t79 * t120 - t78 * t130) * t75, (t79 * t121 - t81 * t130) * t75, 0, 0, 0, 0, 0, (-(-t77 * t138 + t141) * t81 + t77 * t137) * t135 + t3 * qJD(3) + t8 * qJD(4), ((t80 * t138 + t144) * t81 + t80 * t137) * t135 + t4 * qJD(3) + t7 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78 * t111 - t125, t41 * qJD(3) - t81 * t111, 0, 0, 0, 0, 0, t3 * qJD(2) + t12 * qJD(4) - t80 * t125, t4 * qJD(2) + t14 * qJD(4) + t77 * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 * qJD(2) + t12 * qJD(3) + t92 * qJD(4), t7 * qJD(2) + t14 * qJD(3) + t93 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * qJD(4), t6 * qJD(4); 0, 0, 0, 0, t108, t64 * qJD(3), 0, 0, 0, -pkin(2) * t133, -pkin(2) * t131, t73 * t108 - t72 * t110, -t51 * qJD(4) - 0.2e1 * t81 * t100, -t53 * qJD(3) + t78 * t113, t52 * qJD(3) + t78 * t112, -t108, -t9 * qJD(3) - t22 * qJD(4), t10 * qJD(3) + t21 * qJD(4); 0, 0, 0, 0, t107, t122, t131, -t133, 0, -pkin(6) * t131 - t116, pkin(6) * t133 - t115, -t40 + (t73 * t121 + t109) * t81, -0.2e1 * t78 * t110 + t153 * t81, t77 * t133 - t123, t106 + t124, -t46, -t136 + (t98 * t77 - t118) * qJD(3) + t16 * qJD(4), t126 + (t98 * t80 + t119) * qJD(3) + t15 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t85, t77 * t95, t89, t69, t16 * qJD(3) - t35 * qJD(4) - t96, t15 * qJD(3) + t34 * qJD(4) - t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t107, -t122, 0, 0, 0, t116, t115, -t73 * t107 - t40, 0.2e1 * t77 * t89, -t112 + t123, t113 - t124, t46, t24 * qJD(4) + t136, -t23 * qJD(4) - t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t63 * qJD(4), 0, 0, 0, -pkin(3) * t129, -pkin(3) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t153, -t99 * t80, t99 * t77, -t121 / 0.2e1, -pkin(7) * t128 - t90, pkin(7) * t129 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * qJD(2), -t6 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t85, (-t77 * t121 + t132) * t81, (-t114 - t134) * t81, t69, -t24 * qJD(3) + t96, t23 * qJD(3) + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t153, t80 * t120, -t77 * t120, t121 / 0.2e1, t90, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
