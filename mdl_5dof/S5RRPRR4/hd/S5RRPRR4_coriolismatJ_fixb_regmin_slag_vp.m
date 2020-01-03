% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:02:08
% EndTime: 2020-01-03 12:02:11
% DurationCPUTime: 0.92s
% Computational Cost: add. (1197->141), mult. (2448->186), div. (0->0), fcn. (2352->8), ass. (0->126)
t118 = qJD(4) + qJD(5);
t87 = sin(pkin(9));
t91 = sin(qJ(2));
t81 = t87 * t91 * pkin(1);
t94 = cos(qJ(2));
t160 = t94 * pkin(1);
t85 = pkin(2) + t160;
t88 = cos(pkin(9));
t107 = t88 * t85 - t81;
t63 = -pkin(3) - t107;
t83 = -t88 * pkin(2) - pkin(3);
t175 = t83 / 0.2e1 + t63 / 0.2e1;
t92 = cos(qJ(5));
t93 = cos(qJ(4));
t148 = t92 * t93;
t89 = sin(qJ(5));
t90 = sin(qJ(4));
t152 = t89 * t90;
t71 = -t148 + t152;
t173 = t118 * t71;
t119 = qJD(1) + qJD(2);
t149 = t92 * t90;
t151 = t89 * t93;
t73 = t149 + t151;
t29 = t71 ^ 2 - t73 ^ 2;
t172 = t119 * t29;
t79 = -t90 ^ 2 + t93 ^ 2;
t171 = t119 * t79;
t170 = t119 * t93;
t82 = t87 * pkin(2) + pkin(7);
t158 = pkin(8) + t82;
t69 = t158 * t90;
t70 = t158 * t93;
t169 = t118 * (t89 * t69 - t92 * t70);
t168 = t118 * (t92 * t69 + t89 * t70);
t153 = t88 * t91;
t104 = pkin(1) * t153 + t87 * t85;
t64 = pkin(7) + t104;
t159 = pkin(8) + t64;
t43 = t159 * t90;
t44 = t159 * t93;
t167 = t118 * (t89 * t43 - t92 * t44);
t166 = t118 * (t92 * t43 + t89 * t44);
t68 = t88 * t160 - t81;
t165 = -t68 / 0.2e1;
t162 = pkin(4) * t90;
t161 = t93 * pkin(4);
t49 = t63 - t161;
t157 = t49 * t71;
t156 = t49 * t73;
t76 = t83 - t161;
t155 = t76 * t71;
t154 = t76 * t73;
t150 = t90 * t68;
t147 = -t157 / 0.2e1 - t155 / 0.2e1;
t109 = t93 * t165;
t146 = t92 * t109 + t89 * t150 / 0.2e1;
t145 = pkin(1) * qJD(1);
t144 = pkin(1) * qJD(2);
t143 = pkin(4) * qJD(5);
t142 = qJD(4) * pkin(4);
t139 = qJD(1) * t71;
t138 = qJD(1) * t73;
t137 = qJD(1) * t90;
t136 = qJD(1) * t93;
t67 = (t87 * t94 + t153) * pkin(1);
t135 = qJD(2) * t67;
t134 = qJD(2) * t76;
t133 = qJD(2) * t90;
t132 = qJD(2) * t93;
t14 = t104 * t68 - t107 * t67;
t129 = t14 * qJD(1);
t65 = t71 * t162;
t21 = t65 + t156;
t124 = t21 * qJD(1);
t117 = t73 * t162;
t22 = t117 - t157;
t123 = t22 * qJD(1);
t122 = t71 * qJD(5);
t121 = t73 * qJD(5);
t120 = t90 * qJD(4);
t86 = t93 * qJD(4);
t116 = t49 * t139;
t115 = t49 * t138;
t114 = t63 * t137;
t113 = t63 * t136;
t112 = t67 * t139;
t111 = t67 * t138;
t110 = t67 * t137;
t108 = t76 / 0.2e1 + t49 / 0.2e1;
t106 = pkin(1) * t119;
t105 = pkin(4) * t118;
t40 = t118 * t73;
t103 = t165 - t175;
t102 = t117 + t147;
t95 = (-t151 / 0.2e1 - t149 / 0.2e1) * t68;
t5 = -t108 * t73 + t95;
t1 = -t65 + t5;
t27 = t65 + t154;
t101 = t1 * qJD(1) - t27 * qJD(2);
t28 = t117 - t155;
t3 = (t148 / 0.2e1 - t152 / 0.2e1) * t68 + t102;
t100 = -t3 * qJD(1) - t28 * qJD(2);
t99 = t5 * qJD(1) - t73 * t134;
t6 = t108 * t71 + t146;
t98 = t6 * qJD(1) + t71 * t134;
t23 = t103 * t90;
t97 = t23 * qJD(1) - t83 * t133;
t24 = t103 * t93;
t96 = t24 * qJD(1) - t83 * t132;
t8 = t154 / 0.2e1 + t156 / 0.2e1 + t95;
t80 = t90 * t86;
t78 = t79 * qJD(4);
t66 = t90 * t170;
t56 = t67 * t133;
t42 = t73 * t135;
t41 = t71 * t135;
t26 = t93 * t175 + t109;
t25 = -t150 / 0.2e1 + t90 * t175;
t18 = t119 * t73 * t71;
t17 = t71 * t40;
t11 = t118 * t29;
t7 = t146 + t147;
t4 = t102 + t146;
t2 = t65 + t8;
t9 = [0, 0, 0, 0, -t91 * t144, -t94 * t144, t14 * qJD(2), t80, t78, 0, 0, 0, t63 * t120 - t67 * t132, t63 * t86 + t56, -t17, t11, 0, 0, 0, t21 * qJD(4) + t49 * t121 + t41, t22 * qJD(4) - t49 * t122 + t42; 0, 0, 0, 0, -t91 * t106, -t94 * t106, t129 + (-t67 * t88 + t68 * t87) * qJD(2) * pkin(2), t80, t78, 0, 0, 0, t25 * qJD(4) - t67 * t170, t26 * qJD(4) + t110 + t56, -t17, t11, 0, 0, 0, t2 * qJD(4) + t8 * qJD(5) + t112 + t41, t4 * qJD(4) + t7 * qJD(5) + t111 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t66, t171, t86, -t120, 0, t25 * qJD(2) - t64 * t86 + t114, t26 * qJD(2) + t64 * t120 + t113, -t18, t172, -t173, -t40, 0, t2 * qJD(2) + t124 + t167, t4 * qJD(2) + t123 + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t172, -t173, -t40, 0, t8 * qJD(2) + t115 + t167, t7 * qJD(2) - t116 + t166; 0, 0, 0, 0, t91 * t145, t94 * t145, -t129, t80, t78, 0, 0, 0, -t23 * qJD(4) + t67 * t136, -t24 * qJD(4) - t110, -t17, t11, 0, 0, 0, -t1 * qJD(4) - t5 * qJD(5) - t112, t3 * qJD(4) - t6 * qJD(5) - t111; 0, 0, 0, 0, 0, 0, 0, t80, t78, 0, 0, 0, t83 * t120, t83 * t86, -t17, t11, 0, 0, 0, t27 * qJD(4) + t76 * t121, t28 * qJD(4) - t76 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t66, t171, t86, -t120, 0, -t82 * t86 - t97, t82 * t120 - t96, -t18, t172, -t173, -t40, 0, -t101 + t169, -t100 + t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t172, -t173, -t40, 0, -t99 + t169, -t98 + t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t86, 0, 0, 0, 0, 0, -t40, t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t173; 0, 0, 0, 0, 0, 0, 0, -t66, -t171, 0, 0, 0, t23 * qJD(2) - t114, t24 * qJD(2) - t113, t18, -t172, 0, 0, 0, t1 * qJD(2) - t124, -t3 * qJD(2) - t123; 0, 0, 0, 0, 0, 0, 0, -t66, -t171, 0, 0, 0, t97, t96, t18, -t172, 0, 0, 0, t101, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89 * t143, -t92 * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89 * t105, -t92 * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t172, 0, 0, 0, t5 * qJD(2) - t115, t6 * qJD(2) + t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t172, 0, 0, 0, t99, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t142, t92 * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t9;
