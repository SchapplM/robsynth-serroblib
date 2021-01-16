% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:15
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:14:44
% EndTime: 2021-01-15 16:14:46
% DurationCPUTime: 0.70s
% Computational Cost: add. (722->133), mult. (1499->169), div. (0->0), fcn. (1040->4), ass. (0->117)
t95 = cos(qJ(3));
t146 = t95 * pkin(2);
t81 = -pkin(3) - t146;
t152 = t81 / 0.2e1 - pkin(3) / 0.2e1;
t122 = qJD(2) + qJD(3);
t92 = sin(qJ(4));
t90 = t92 ^ 2;
t94 = cos(qJ(4));
t91 = t94 ^ 2;
t77 = t90 + t91;
t154 = t122 * t77;
t78 = t91 - t90;
t153 = t122 * t78;
t149 = pkin(4) * t92;
t93 = sin(qJ(3));
t148 = t93 * pkin(2);
t147 = t94 * pkin(4);
t80 = pkin(7) + t148;
t132 = qJ(5) + t80;
t48 = t132 * t92;
t145 = t48 * t92;
t49 = t132 * t94;
t144 = t49 * t94;
t82 = -pkin(3) - t147;
t65 = t82 - t146;
t56 = t65 * t92;
t143 = t65 * t94;
t139 = qJ(5) + pkin(7);
t66 = t139 * t92;
t142 = t66 * t92;
t67 = t139 * t94;
t141 = t67 * t94;
t71 = t82 * t92;
t140 = t82 * t94;
t47 = t77 * t146;
t68 = t77 * qJD(5);
t138 = t47 * qJD(3) + t68;
t136 = pkin(2) * qJD(2);
t117 = t93 * t136;
t74 = t92 * t117;
t135 = pkin(2) * qJD(3);
t120 = t93 * t135;
t76 = t92 * t120;
t137 = t74 + t76;
t134 = pkin(3) * qJD(3);
t17 = t144 + t145;
t8 = (t17 * t95 + t65 * t93) * pkin(2);
t133 = t8 * qJD(2);
t131 = qJD(2) * t81;
t130 = t17 * qJD(2);
t121 = t92 * t147;
t30 = -t56 + t121;
t129 = t30 * qJD(2);
t89 = t90 * pkin(4);
t40 = t89 + t143;
t128 = t40 * qJD(2);
t127 = t47 * qJD(2);
t126 = t49 * qJD(4);
t125 = t67 * qJD(4);
t87 = t92 * qJD(4);
t124 = t92 * qJD(5);
t88 = t94 * qJD(4);
t123 = t94 * qJD(5);
t119 = pkin(4) * t88;
t118 = pkin(4) * t124;
t115 = t148 / 0.2e1;
t114 = -t146 / 0.2e1;
t113 = t92 * t131;
t112 = t94 * t131;
t111 = -t71 / 0.2e1 - t56 / 0.2e1;
t110 = pkin(2) * t122;
t109 = t94 * t120;
t108 = t92 * t114;
t107 = t94 * t114;
t63 = t122 * t92;
t64 = t122 * t94;
t21 = t141 + t142;
t38 = t144 / 0.2e1;
t13 = t38 - t144 / 0.2e1;
t3 = pkin(4) * t56;
t106 = t13 * qJD(1) + t3 * qJD(2);
t4 = t115 + (-t67 / 0.2e1 - t49 / 0.2e1) * t94 + (-t66 / 0.2e1 - t48 / 0.2e1) * t92;
t105 = -t4 * qJD(2) + t21 * qJD(3);
t75 = t94 * t117;
t104 = -t75 - t109;
t11 = (t114 + t147) * t92 + t111;
t41 = -t71 + t121;
t103 = -t11 * qJD(2) - t41 * qJD(3);
t53 = t141 / 0.2e1;
t19 = t53 - t141 / 0.2e1;
t102 = -t13 * qJD(2) - t19 * qJD(3);
t99 = t114 - t82 / 0.2e1 - t65 / 0.2e1;
t14 = t99 * t94 - t89;
t55 = t89 + t140;
t101 = -t14 * qJD(2) + t55 * qJD(3);
t100 = t114 - t152;
t26 = t100 * t92;
t98 = t26 * qJD(2) + t92 * t134;
t27 = t100 * t94;
t97 = t27 * qJD(2) + t94 * t134;
t1 = t99 * t149;
t9 = pkin(4) * t71;
t96 = t19 * qJD(1) - t1 * qJD(2) + t9 * qJD(3);
t86 = pkin(4) * t87;
t79 = t92 * t88;
t69 = t78 * qJD(4);
t57 = pkin(4) * t63;
t44 = t92 * t64;
t29 = t152 * t94 + t107;
t28 = t152 * t92 + t108;
t18 = t19 * qJD(4);
t15 = t89 + t140 / 0.2e1 + t143 / 0.2e1 + t107;
t12 = (t114 - t147) * t92 - t111;
t10 = t13 * qJD(4);
t5 = t53 + t38 + t142 / 0.2e1 + t145 / 0.2e1 + t115;
t2 = pkin(4) * t108 + (t65 + t82) * t149 / 0.2e1;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t88, -t87, -t88, 0, -t102 - t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, -t120, -t95 * t135, t79, t69, 0, 0, 0, t81 * t87 - t109, t81 * t88 + t76, -t30 * qJD(4) - t109, t40 * qJD(4) + t76, t138, t8 * qJD(3) + t3 * qJD(4) + t17 * qJD(5); 0, 0, 0, 0, 0, -t93 * t110, -t95 * t110, t79, t69, 0, 0, 0, t28 * qJD(4) + t104, t29 * qJD(4) + t137, t12 * qJD(4) + t104, t15 * qJD(4) + t137, t127 + t138, t133 + t2 * qJD(4) + t5 * qJD(5) + (t21 * t95 + t82 * t93) * t135; 0, 0, 0, 0, 0, 0, 0, t44, t153, t88, -t87, 0, t28 * qJD(3) - t80 * t88 + t113, t29 * qJD(3) + t80 * t87 + t112, t12 * qJD(3) - t126 - t129, t15 * qJD(3) + t48 * qJD(4) + t128, -t119, -pkin(4) * t126 + t2 * qJD(3) + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t5 * qJD(3) + t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, t117, t95 * t136, t79, t69, 0, 0, 0, -t26 * qJD(4) + t75, -t27 * qJD(4) - t74, -t11 * qJD(4) + t75, -t14 * qJD(4) - t74, t68 - t127, -t1 * qJD(4) - t4 * qJD(5) - t133; 0, 0, 0, 0, 0, 0, 0, t79, t69, 0, 0, 0, -pkin(3) * t87, -pkin(3) * t88, -t41 * qJD(4), t55 * qJD(4), t68, t9 * qJD(4) + t21 * qJD(5); 0, 0, 0, 0, 0, 0, 0, t44, t153, t88, -t87, 0, -pkin(7) * t88 - t98, pkin(7) * t87 - t97, t103 - t125, t66 * qJD(4) + t101, -t119, -pkin(4) * t125 + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102; 0, 0, 0, 0, 0, 0, 0, -t44, -t153, 0, 0, 0, t26 * qJD(3) - t113, t27 * qJD(3) - t112, t11 * qJD(3) - t124 + t129, t14 * qJD(3) - t123 - t128, 0, t1 * qJD(3) - t106 - t118; 0, 0, 0, 0, 0, 0, 0, -t44, -t153, 0, 0, 0, t98, t97, -t103 - t124, -t101 - t123, 0, -t96 - t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t64, 0, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t88, -t154, t4 * qJD(3) - t130 + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t88, -t154, -t105 + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t64, 0, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
