% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:50
% EndTime: 2019-12-31 18:02:54
% DurationCPUTime: 1.11s
% Computational Cost: add. (536->155), mult. (1228->248), div. (0->0), fcn. (1065->6), ass. (0->137)
t74 = sin(qJ(5));
t168 = 0.2e1 * t74;
t77 = cos(qJ(4));
t123 = t77 * qJD(4);
t75 = sin(qJ(4));
t125 = t75 * qJD(5);
t62 = t74 * t125;
t76 = cos(qJ(5));
t167 = -t76 * t123 + t62;
t127 = t75 * qJD(1);
t117 = t76 * t127;
t68 = t74 ^ 2;
t70 = t76 ^ 2;
t58 = t70 - t68;
t81 = t58 * qJD(4) + t117 * t168;
t73 = cos(pkin(8));
t153 = t76 * t73;
t157 = t74 * t77;
t72 = sin(pkin(8));
t32 = t72 * t157 + t153;
t166 = t32 / 0.2e1;
t152 = t76 * t77;
t159 = t74 * t73;
t33 = t72 * t152 - t159;
t165 = t33 / 0.2e1;
t69 = t75 ^ 2;
t164 = t69 / 0.2e1;
t163 = t75 * pkin(4);
t162 = t77 * pkin(7);
t51 = t162 - t163;
t161 = t74 * t51;
t160 = t74 * t72;
t158 = t74 * t75;
t156 = t76 * t51;
t155 = t76 * t69;
t154 = t76 * t72;
t78 = -pkin(1) - pkin(2);
t151 = t73 * qJ(2) + t72 * t78;
t71 = t77 ^ 2;
t59 = t71 - t69;
t41 = -pkin(6) + t151;
t121 = t41 * t157;
t98 = -t72 * qJ(2) + t73 * t78;
t40 = pkin(3) - t98;
t94 = t77 * pkin(4) + t75 * pkin(7);
t79 = t40 + t94;
t12 = -t76 * t79 + t121;
t30 = t41 * t158;
t1 = (t30 + t156) * t77 + (t12 - 0.2e1 * t121) * t75;
t150 = t1 * qJD(1);
t120 = t75 * t41 * t76;
t119 = t41 * t152;
t13 = t74 * t79 + t119;
t2 = -t13 * t75 + (t120 + t161) * t77;
t149 = t2 * qJD(1);
t3 = -t69 * t41 * t74 - t12 * t77;
t148 = t3 * qJD(1);
t4 = -t13 * t77 - t41 * t155;
t147 = t4 * qJD(1);
t102 = -t159 / 0.2e1;
t118 = 0.1e1 / 0.2e1 + t164;
t7 = (t102 + t165) * t77 + t118 * t154;
t146 = t7 * qJD(1);
t100 = -t153 / 0.2e1;
t8 = (t100 - t32 / 0.2e1) * t77 - t118 * t160;
t145 = t8 * qJD(1);
t144 = qJD(4) * t74;
t143 = qJD(4) * t76;
t142 = qJD(5) * t74;
t141 = qJD(5) * t76;
t140 = qJD(5) * t77;
t101 = -t157 / 0.2e1;
t92 = t72 * t101 + t166;
t14 = (t153 / 0.2e1 + t92) * t75;
t139 = t14 * qJD(1);
t99 = -t152 / 0.2e1;
t91 = t72 * t99 + t165;
t17 = (t102 + t91) * t75;
t138 = t17 * qJD(1);
t18 = t151 * t73 - t98 * t72;
t137 = t18 * qJD(1);
t21 = (-t73 * t157 + t154) * t77 - t69 * t159;
t136 = t21 * qJD(1);
t22 = (t73 * t152 + t160) * t77 + t69 * t153;
t135 = t22 * qJD(1);
t43 = t59 * t74;
t134 = t43 * qJD(1);
t44 = t76 * t71 - t155;
t133 = t44 * qJD(1);
t132 = t59 * qJD(1);
t131 = t72 * qJD(1);
t130 = t72 * qJD(2);
t129 = t73 * qJD(1);
t128 = t73 * qJD(2);
t126 = t75 * qJD(4);
t124 = t77 * qJD(1);
t122 = qJ(2) * qJD(1);
t116 = t74 * t140;
t115 = t76 * t140;
t114 = t73 * t124;
t113 = t74 * t126;
t112 = t74 * t141;
t111 = t74 * t143;
t110 = t75 * t123;
t109 = t72 * t127;
t108 = t72 * t126;
t107 = t73 * t127;
t106 = t75 * t124;
t105 = t76 * t126;
t103 = t72 * t124;
t96 = t74 * t105;
t95 = qJD(5) + t124;
t93 = t95 * t76;
t90 = qJD(1) * t40 - t128;
t89 = t162 / 0.2e1 - t163 / 0.2e1;
t84 = t51 / 0.2e1 + t89;
t19 = t84 * t74;
t88 = pkin(4) * t143 - t19 * qJD(1);
t20 = t84 * t76;
t87 = pkin(4) * t144 + t20 * qJD(1);
t86 = t75 * t93;
t34 = (t68 / 0.2e1 - t70 / 0.2e1) * t75;
t85 = t34 * qJD(1) + t111;
t83 = t74 * qJD(1) * t155 - t34 * qJD(4);
t42 = t58 * t69;
t82 = -t42 * qJD(1) + 0.2e1 * t96;
t80 = t74 * t123 + t76 * t125;
t63 = -t126 / 0.2e1;
t37 = (t124 + qJD(5) / 0.2e1) * t75;
t31 = t34 * qJD(5);
t16 = t73 * t158 / 0.2e1 + t91 * t75;
t15 = (t100 + t92) * t75;
t10 = -t33 * t77 / 0.2e1 + t73 * t101 + (-t69 / 0.2e1 + 0.1e1 / 0.2e1) * t154;
t9 = t77 * t166 + t73 * t99 + (t164 - 0.1e1 / 0.2e1) * t160;
t6 = t30 + t156 / 0.2e1 - t89 * t76;
t5 = t120 - t161 / 0.2e1 + t89 * t74;
t11 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), t130, t128, t18 * qJD(2), t110, t59 * qJD(4), 0, 0, 0, -t40 * t126 + t77 * t130, -t40 * t123 - t75 * t130, t70 * t110 - t69 * t112, -t42 * qJD(5) - 0.2e1 * t77 * t96, -t44 * qJD(4) + t75 * t116, t43 * qJD(4) + t75 * t115, -t110, t21 * qJD(2) + t1 * qJD(4) + t4 * qJD(5), -t22 * qJD(2) - t2 * qJD(4) - t3 * qJD(5); 0, 0, 0, 0, qJD(1), t122, t131, t129, t137, 0, 0, 0, 0, 0, t103, -t109, 0, 0, 0, 0, 0, t15 * qJD(4) + t10 * qJD(5) + t136, t16 * qJD(4) + t9 * qJD(5) - t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t132, -t123, t126, 0, -t41 * t123 - t40 * t127, -t40 * t124 + t41 * t126, t31 + (t70 * t127 - t111) * t77, 0.2e1 * t75 * t112 - t81 * t77, -t113 - t133, -t105 + t134, -t37, t150 + t15 * qJD(2) + (t94 * t74 - t119) * qJD(4) + t6 * qJD(5), -t149 + t16 * qJD(2) + (t94 * t76 + t121) * qJD(4) + t5 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t82, t74 * t106 + t62, t86, t63, t10 * qJD(2) + t6 * qJD(4) - t13 * qJD(5) + t147, t9 * qJD(2) + t5 * qJD(4) + t12 * qJD(5) - t148; 0, 0, 0, 0, -qJD(1), -t122, -t131, -t129, -t137, 0, 0, 0, 0, 0, t73 * t126 - t103, t73 * t123 + t109, 0, 0, 0, 0, 0, t14 * qJD(4) - t7 * qJD(5) - t136, t17 * qJD(4) - t8 * qJD(5) + t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72 * t123 + t107, t108 + t114, 0, 0, 0, 0, 0, t167 * t72 + t139, t80 * t72 + t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33 * qJD(5) + t74 * t108 - t146, t32 * qJD(5) + t72 * t105 - t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t123, 0, 0, 0, 0, 0, -t105 - t116, t113 - t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t132, 0, 0, 0, t90 * t75, t90 * t77, -t70 * t106 + t31, t86 * t168, t115 + t133, -t116 - t134, t37, -t14 * qJD(2) - t20 * qJD(5) - t150, -t17 * qJD(2) + t19 * qJD(5) + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t114, 0, 0, 0, 0, 0, -t139, -t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t58 * qJD(5), 0, 0, 0, -pkin(4) * t142, -pkin(4) * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t81, t93, -t95 * t74, t127 / 0.2e1, -pkin(7) * t141 - t87, pkin(7) * t142 - t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t82, (-t74 * t127 - t143) * t77, (-t117 + t144) * t77, t63, t7 * qJD(2) + t20 * qJD(4) - t147, t8 * qJD(2) - t19 * qJD(4) + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t81, -t76 * t124, t74 * t124, -t127 / 0.2e1, t87, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t11;
