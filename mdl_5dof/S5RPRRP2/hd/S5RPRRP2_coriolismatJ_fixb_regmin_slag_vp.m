% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:36:31
% EndTime: 2021-01-15 12:36:34
% DurationCPUTime: 0.74s
% Computational Cost: add. (1166->138), mult. (2174->173), div. (0->0), fcn. (1709->6), ass. (0->123)
t125 = qJD(1) + qJD(3);
t98 = sin(qJ(4));
t95 = t98 ^ 2;
t100 = cos(qJ(4));
t96 = t100 ^ 2;
t86 = t95 + t96;
t159 = t125 * t86;
t87 = t96 - t95;
t158 = t125 * t87;
t157 = -pkin(3) / 0.2e1;
t112 = cos(pkin(8)) * pkin(1) + pkin(2);
t150 = cos(qJ(3));
t154 = pkin(1) * sin(pkin(8));
t99 = sin(qJ(3));
t62 = -t150 * t112 + t154 * t99;
t156 = -t62 / 0.2e1;
t155 = t98 / 0.2e1;
t153 = pkin(4) * t98;
t152 = t95 * pkin(4);
t151 = t100 / 0.2e1;
t149 = t100 * pkin(4);
t105 = t99 * t112;
t114 = t150 * t154;
t63 = t114 + t105;
t61 = pkin(7) + t63;
t140 = qJ(5) + t61;
t40 = t140 * t98;
t148 = t40 * t98;
t60 = -pkin(3) + t62;
t48 = t60 - t149;
t44 = t48 * t98;
t146 = pkin(7) + qJ(5);
t79 = t146 * t98;
t147 = t79 * t98;
t89 = -pkin(3) - t149;
t84 = t89 * t98;
t23 = t86 * t62;
t82 = t86 * qJD(5);
t145 = -t23 * qJD(3) + t82;
t130 = t98 * qJD(3);
t51 = t63 * t130;
t131 = t98 * qJD(1);
t52 = t63 * t131;
t144 = t51 + t52;
t41 = t140 * t100;
t142 = t41 * t100;
t16 = t142 + t148;
t3 = -t16 * t62 + t48 * t63;
t143 = t3 * qJD(1);
t80 = t146 * t100;
t141 = t80 * t100;
t139 = t16 * qJD(1);
t138 = t23 * qJD(1);
t124 = t98 * t149;
t28 = -t44 + t124;
t137 = t28 * qJD(1);
t34 = t48 * t100 + t152;
t136 = t34 * qJD(1);
t135 = t41 * qJD(4);
t134 = t62 * qJD(1);
t133 = t63 * qJD(1);
t132 = t80 * qJD(4);
t93 = t98 * qJD(4);
t129 = t98 * qJD(5);
t128 = t100 * qJD(1);
t127 = t100 * qJD(3);
t94 = t100 * qJD(4);
t126 = t100 * qJD(5);
t123 = pkin(4) * t129;
t121 = pkin(3) / 0.2e1 - t60 / 0.2e1;
t120 = pkin(4) * t94;
t119 = t60 * t131;
t49 = t62 * t155;
t118 = t84 / 0.2e1 + t44 / 0.2e1;
t117 = t89 / 0.2e1 + t48 / 0.2e1;
t116 = t60 * t128;
t115 = t63 * t127;
t77 = t125 * t98;
t78 = t125 * t100;
t113 = t156 + t117;
t53 = t63 * t128;
t111 = -t53 - t115;
t43 = t141 + t147;
t31 = t142 / 0.2e1;
t14 = t31 - t142 / 0.2e1;
t4 = pkin(4) * t44;
t110 = t4 * qJD(1) + t14 * qJD(2);
t101 = t114 / 0.2e1 + t105 / 0.2e1;
t5 = (-t79 / 0.2e1 - t40 / 0.2e1) * t98 + (-t80 / 0.2e1 - t41 / 0.2e1) * t100 + t101;
t109 = -t5 * qJD(1) + t43 * qJD(3);
t64 = -t84 + t124;
t9 = (-t149 + t156) * t98 + t118;
t108 = t9 * qJD(1) - t64 * qJD(3);
t12 = t100 * t113 + t152;
t71 = t89 * t100 + t152;
t107 = t12 * qJD(1) + t71 * qJD(3);
t69 = t141 / 0.2e1;
t37 = t69 - t141 / 0.2e1;
t106 = -t14 * qJD(1) - t37 * qJD(3);
t18 = t121 * t98 + t49;
t104 = pkin(3) * t130 + t18 * qJD(1);
t50 = t62 * t151;
t19 = t100 * t121 + t50;
t103 = pkin(3) * t127 + t19 * qJD(1);
t1 = t113 * t153;
t17 = pkin(4) * t84;
t102 = t1 * qJD(1) + t37 * qJD(2) + t17 * qJD(3);
t92 = pkin(4) * t93;
t88 = t98 * t94;
t83 = t87 * qJD(4);
t72 = pkin(4) * t77;
t65 = t98 * t78;
t59 = t63 * qJD(3);
t58 = t62 * qJD(3);
t35 = t37 * qJD(4);
t21 = t100 * t157 + t60 * t151 + t50;
t20 = t60 * t155 + t157 * t98 + t49;
t13 = t100 * t117 + t152 + t50;
t11 = t14 * qJD(4);
t10 = t118 + t49 - t124;
t6 = t69 + t31 + t147 / 0.2e1 + t148 / 0.2e1 + t101;
t2 = pkin(4) * t49 + (t48 + t89) * t153 / 0.2e1;
t7 = [0, 0, 0, 0, 0, -t59, t58, t88, t83, 0, 0, 0, t60 * t93 - t115, t60 * t94 + t51, -t28 * qJD(4) - t115, t34 * qJD(4) + t51, t145, t3 * qJD(3) + t4 * qJD(4) + t16 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, -t59 - t133, t58 + t134, t88, t83, 0, 0, 0, t20 * qJD(4) + t111, t21 * qJD(4) + t144, t10 * qJD(4) + t111, t13 * qJD(4) + t144, -t138 + t145, t143 + (-t43 * t62 + t63 * t89) * qJD(3) + t2 * qJD(4) + t6 * qJD(5); 0, 0, 0, 0, 0, 0, 0, t65, t158, t94, -t93, 0, t20 * qJD(3) - t61 * t94 + t119, t21 * qJD(3) + t61 * t93 + t116, t10 * qJD(3) - t135 - t137, t13 * qJD(3) + t40 * qJD(4) + t136, -t120, -pkin(4) * t135 + t2 * qJD(3) + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t6 * qJD(3) + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, -t94, -t93, -t94, 0, -t106 - t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t133, -t134, t88, t83, 0, 0, 0, -t18 * qJD(4) + t53, -t19 * qJD(4) - t52, t9 * qJD(4) + t53, t12 * qJD(4) - t52, t82 + t138, t1 * qJD(4) - t5 * qJD(5) - t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, t88, t83, 0, 0, 0, -pkin(3) * t93, -pkin(3) * t94, -t64 * qJD(4), t71 * qJD(4), t82, t17 * qJD(4) + t43 * qJD(5); 0, 0, 0, 0, 0, 0, 0, t65, t158, t94, -t93, 0, -pkin(7) * t94 - t104, pkin(7) * t93 - t103, t108 - t132, t79 * qJD(4) + t107, -t120, -pkin(4) * t132 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t109; 0, 0, 0, 0, 0, 0, 0, -t65, -t158, 0, 0, 0, t18 * qJD(3) - t119, t19 * qJD(3) - t116, -t9 * qJD(3) - t129 + t137, -t12 * qJD(3) - t126 - t136, 0, -t1 * qJD(3) - t110 - t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106; 0, 0, 0, 0, 0, 0, 0, -t65, -t158, 0, 0, 0, t104, t103, -t108 - t129, -t107 - t126, 0, -t102 - t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t78, 0, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t94, -t159, t5 * qJD(3) - t139 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t94, -t159, -t109 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t78, 0, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
