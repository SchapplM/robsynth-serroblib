% Calculate minimal parameter regressor of coriolis matrix for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x19]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:45
% EndTime: 2019-12-31 17:06:48
% DurationCPUTime: 0.85s
% Computational Cost: add. (999->122), mult. (2197->213), div. (0->0), fcn. (2281->6), ass. (0->122)
t83 = sin(qJ(4));
t165 = 0.2e1 * t83;
t142 = cos(pkin(7));
t84 = sin(qJ(2));
t105 = t142 * t84;
t150 = -qJ(3) - pkin(5);
t86 = cos(qJ(2));
t71 = t150 * t86;
t82 = sin(pkin(7));
t43 = -t150 * t105 - t82 * t71;
t164 = t43 / 0.2e1;
t154 = t82 * t84;
t87 = -t142 * t71 + t150 * t154;
t158 = t87 * t83;
t85 = cos(qJ(4));
t157 = t87 * t85;
t65 = -t142 * t86 + t154;
t63 = t65 ^ 2;
t153 = t82 * t86;
t67 = t105 + t153;
t64 = t67 ^ 2;
t163 = -t64 - t63;
t31 = t83 * t67;
t106 = 0.2e1 * t85 * t31;
t80 = t83 ^ 2;
t81 = t85 ^ 2;
t74 = t81 - t80;
t91 = qJD(1) * t106 - t74 * qJD(2);
t160 = -t67 / 0.2e1;
t159 = t84 * pkin(2);
t156 = t43 * t85;
t155 = t65 * t67;
t34 = t67 * pkin(3) + t65 * pkin(6) + t159;
t152 = t83 * t34;
t151 = t85 * t34;
t148 = qJD(2) * pkin(2);
t107 = -t86 * pkin(2) - pkin(1);
t89 = t65 * pkin(3) - t67 * pkin(6) + t107;
t12 = -t85 * t89 + t158;
t1 = (-t12 + t158) * t67 + t151 * t65;
t147 = t1 * qJD(1);
t13 = t83 * t89 + t157;
t2 = (-t13 + t157) * t67 - t152 * t65;
t146 = t2 * qJD(1);
t8 = t12 * t65 - t43 * t31;
t144 = t8 * qJD(1);
t9 = -t13 * t65 + t67 * t156;
t143 = t9 * qJD(1);
t141 = qJD(1) * t67;
t140 = qJD(1) * t85;
t139 = qJD(1) * t86;
t138 = qJD(2) * t83;
t137 = qJD(2) * t84;
t136 = qJD(2) * t85;
t135 = qJD(2) * t86;
t134 = qJD(3) * t85;
t133 = qJD(4) * t83;
t132 = qJD(4) * t85;
t10 = t107 * t159;
t131 = t10 * qJD(1);
t11 = t43 * t67 - t65 * t87;
t130 = t11 * qJD(1);
t118 = t64 - t63;
t14 = t118 * t83;
t129 = t14 * qJD(1);
t15 = t163 * t83;
t128 = t15 * qJD(1);
t16 = t118 * t85;
t127 = t16 * qJD(1);
t88 = -t82 * t65 / 0.2e1 + t142 * t160;
t19 = (-t84 / 0.2e1 + t88) * pkin(2);
t126 = t19 * qJD(1);
t28 = t83 * t65;
t125 = t28 * qJD(1);
t124 = t31 * qJD(1);
t33 = t85 * t65;
t123 = t33 * qJD(1);
t36 = t163 * t85;
t122 = t36 * qJD(1);
t121 = t163 * qJD(1);
t62 = t105 / 0.2e1 + t153 / 0.2e1;
t120 = t62 * qJD(1);
t75 = -t84 ^ 2 + t86 ^ 2;
t119 = t75 * qJD(1);
t117 = pkin(1) * t84 * qJD(1);
t116 = pkin(1) * t139;
t115 = t65 * t132;
t114 = t65 * t141;
t113 = qJD(2) * t155;
t112 = t83 * t132;
t111 = t83 * t136;
t110 = t84 * t139;
t109 = t67 * t140;
t108 = t164 - t43 / 0.2e1;
t104 = -qJD(1) * t65 - qJD(4);
t102 = qJD(2) * t106;
t77 = t82 * pkin(2) + pkin(6);
t78 = -t142 * pkin(2) - pkin(3);
t101 = -t65 * t78 - t67 * t77;
t100 = t104 * t85;
t99 = t77 * t65 / 0.2e1 + t78 * t160;
t90 = t34 / 0.2e1 + t99;
t4 = t108 * t85 + t90 * t83;
t98 = -t4 * qJD(1) - t78 * t136;
t6 = t108 * t83 - t90 * t85;
t97 = -t6 * qJD(1) - t78 * t138;
t96 = t67 * t100;
t27 = (t80 / 0.2e1 - t81 / 0.2e1) * t67;
t95 = -t27 * qJD(1) + t111;
t94 = t62 * qJD(4) + t114;
t93 = t64 * t83 * t140 + t27 * qJD(2);
t35 = t74 * t64;
t92 = t35 * qJD(1) + t102;
t61 = t62 * qJD(2);
t60 = t67 * t136;
t24 = t28 * qJD(4);
t23 = t27 * qJD(4);
t18 = t159 / 0.2e1 + t88 * pkin(2);
t17 = -t125 - t133;
t7 = t151 / 0.2e1 - t99 * t85 + t164 * t165;
t5 = t156 / 0.2e1 + t85 * t164 - t152 / 0.2e1 + t99 * t83;
t3 = [0, 0, 0, t84 * t135, t75 * qJD(2), 0, 0, 0, -pkin(1) * t137, -pkin(1) * t135, -t163 * qJD(3), t10 * qJD(2) + t11 * qJD(3), -t64 * t112 - t81 * t113, -t35 * qJD(4) + t102 * t65, t16 * qJD(2) - t133 * t155, -t14 * qJD(2) - t115 * t67, t113, t1 * qJD(2) - t15 * qJD(3) + t9 * qJD(4), t2 * qJD(2) - t36 * qJD(3) + t8 * qJD(4); 0, 0, 0, t110, t119, t135, -t137, 0, -pkin(5) * t135 - t117, pkin(5) * t137 - t116, (t142 * t65 - t67 * t82) * t148, t131 + (-t142 * t87 - t43 * t82) * t148 + t18 * qJD(3), -t23 + (-t81 * t141 - t111) * t65, -0.2e1 * t67 * t112 + t91 * t65, t138 * t67 + t127, t60 - t129, t94, t147 + (t101 * t83 - t157) * qJD(2) + t7 * qJD(4), t146 + (t101 * t85 + t158) * qJD(2) + t5 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t18 * qJD(2) + t130, 0, 0, 0, 0, 0, -t128, -t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t92, t104 * t31, t96, t61, t7 * qJD(2) - t13 * qJD(4) + t143, t5 * qJD(2) + t12 * qJD(4) + t144; 0, 0, 0, -t110, -t119, 0, 0, 0, t117, t116, 0, t19 * qJD(3) - t131, t81 * t114 - t23, t96 * t165, t33 * qJD(4) - t127, -t24 + t129, -t94, t6 * qJD(4) - t134 * t67 - t147, t31 * qJD(3) + t4 * qJD(4) - t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t74 * qJD(4), 0, 0, 0, t78 * t133, t78 * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, 0, 0, 0, 0, 0, -t109, t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t91, t123 + t132, t17, -t120, -t132 * t77 - t97, t133 * t77 - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t121, -t19 * qJD(2) - t130, 0, 0, 0, 0, 0, -t24 + t60 + t128, -t31 * qJD(2) - t115 + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, 0, 0, 0, 0, 0, t109, -t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t92, -t33 * qJD(2) + t114 * t83, t28 * qJD(2) + t109 * t65, t61, -t6 * qJD(2) + t28 * qJD(3) - t143, -t4 * qJD(2) + t134 * t65 - t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t91, -t123, t125, t120, t97, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t65 * t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
