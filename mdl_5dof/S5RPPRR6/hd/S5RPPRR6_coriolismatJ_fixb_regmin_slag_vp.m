% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:10
% EndTime: 2019-12-31 17:58:12
% DurationCPUTime: 0.87s
% Computational Cost: add. (1117->114), mult. (2398->177), div. (0->0), fcn. (2569->8), ass. (0->112)
t88 = sin(pkin(9));
t92 = sin(qJ(4));
t145 = t92 * t88;
t149 = cos(qJ(4));
t89 = cos(pkin(9));
t74 = -t149 * t89 + t145;
t91 = sin(qJ(5));
t36 = t91 * t74;
t29 = t36 * qJD(5);
t113 = t149 * t88;
t144 = t92 * t89;
t76 = t113 + t144;
t70 = t76 * qJD(4);
t93 = cos(qJ(5));
t66 = t93 * t70;
t155 = -t66 + t29;
t72 = t74 ^ 2;
t73 = t76 ^ 2;
t154 = t72 + t73;
t121 = t73 - t72;
t146 = t91 * t93;
t112 = 0.2e1 * t76 * t146;
t86 = t91 ^ 2;
t87 = t93 ^ 2;
t81 = t87 - t86;
t95 = qJD(1) * t112 - t81 * qJD(4);
t151 = t76 * pkin(4);
t152 = t74 * pkin(7);
t49 = t151 + t152;
t153 = t49 / 0.2e1;
t82 = sin(pkin(8)) * pkin(1) + qJ(3);
t150 = pkin(6) + t82;
t68 = t150 * t89;
t24 = t150 * t113 + t92 * t68;
t148 = t24 * t93;
t25 = -t150 * t145 + t149 * t68;
t147 = t91 * t25;
t37 = t91 * t76;
t143 = t93 * t25;
t40 = t93 * t74;
t78 = t88 ^ 2 + t89 ^ 2;
t107 = t74 * pkin(4) - t76 * pkin(7);
t77 = -cos(pkin(8)) * pkin(1) - t89 * pkin(3) - pkin(2);
t94 = t107 + t77;
t11 = -t93 * t94 + t147;
t1 = t49 * t40 + (-t11 + t147) * t76;
t142 = t1 * qJD(1);
t8 = t11 * t74 - t24 * t37;
t141 = t8 * qJD(1);
t12 = t91 * t94 + t143;
t9 = -t12 * t74 + t76 * t148;
t140 = t9 * qJD(1);
t139 = qJD(3) * t93;
t138 = qJD(4) * t93;
t137 = qJD(5) * t91;
t136 = qJD(5) * t93;
t14 = t121 * t91;
t134 = t14 * qJD(1);
t15 = t154 * t91;
t133 = t15 * qJD(1);
t16 = t121 * t93;
t132 = t16 * qJD(1);
t131 = t121 * qJD(1);
t28 = t36 * qJD(1);
t130 = t37 * qJD(1);
t129 = t40 * qJD(1);
t128 = t40 * qJD(4);
t44 = t154 * t93;
t127 = t44 * qJD(1);
t50 = t78 * t82;
t126 = t50 * qJD(1);
t71 = t113 / 0.2e1 + t144 / 0.2e1;
t125 = t71 * qJD(1);
t124 = t74 * qJD(1);
t69 = t74 * qJD(4);
t123 = t76 * qJD(1);
t122 = t78 * qJD(1);
t120 = t76 * t137;
t119 = t76 * t136;
t118 = t74 * t123;
t117 = t74 * t70;
t116 = t91 * t136;
t115 = t91 * t138;
t114 = t93 * t123;
t111 = qJD(1) * t77 + qJD(3);
t110 = -qJD(5) - t124;
t108 = qJD(4) * t112;
t2 = -t49 * t36 + (-t12 + t143) * t76;
t106 = t2 * qJD(1);
t105 = t110 * t93;
t104 = t152 / 0.2e1 + t151 / 0.2e1;
t98 = t153 + t104;
t4 = t98 * t91;
t103 = pkin(4) * t138 - t4 * qJD(1);
t6 = t98 * t93;
t102 = pkin(4) * t91 * qJD(4) + t6 * qJD(1);
t101 = t76 * t105;
t34 = (t86 / 0.2e1 - t87 / 0.2e1) * t76;
t100 = -t34 * qJD(1) + t115;
t99 = t71 * qJD(5) + t118;
t97 = t73 * qJD(1) * t146 + t34 * qJD(4);
t43 = t81 * t73;
t96 = t43 * qJD(1) + t108;
t67 = t71 * qJD(4);
t65 = t91 * t70;
t33 = t40 * qJD(5);
t27 = t36 * qJD(4);
t26 = t34 * qJD(5);
t18 = -t28 - t137;
t7 = t24 * t91 + (-t104 + t153) * t93;
t5 = t148 + (-t49 / 0.2e1 + t104) * t91;
t3 = [0, 0, 0, 0, 0, 0, t78 * qJD(3), t50 * qJD(3), -t117, -t121 * qJD(4), 0, 0, 0, t77 * t70, -t77 * t69, -t73 * t116 - t87 * t117, -t43 * qJD(5) + t74 * t108, t16 * qJD(4) - t74 * t120, -t14 * qJD(4) - t74 * t119, t117, t15 * qJD(3) + t1 * qJD(4) + t9 * qJD(5), t44 * qJD(3) + t2 * qJD(4) + t8 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t122, t126, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, t127; 0, 0, 0, 0, 0, 0, 0, 0, -t118, -t131, -t69, -t70, 0, -t25 * qJD(4) + t77 * t123, t24 * qJD(4) - t77 * t124, -t26 + (-t87 * t123 - t115) * t74, -0.2e1 * t76 * t116 + t95 * t74, t65 + t132, t66 - t134, t99, t142 + (t107 * t91 - t143) * qJD(4) + t7 * qJD(5), (t107 * t93 + t147) * qJD(4) + t5 * qJD(5) + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96, t110 * t37, t101, t67, t7 * qJD(4) - t12 * qJD(5) + t140, t5 * qJD(4) + t11 * qJD(5) + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t69, 0, 0, 0, 0, 0, t155, t33 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 - t119, t120 + t128; 0, 0, 0, 0, 0, 0, -t122, -t126, 0, 0, 0, 0, 0, t70, -t69, 0, 0, 0, 0, 0, -t133 - t155, -t37 * qJD(4) - t74 * t136 - t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, -t124, 0, 0, 0, 0, 0, t114, -t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t105; 0, 0, 0, 0, 0, 0, 0, 0, t118, t131, 0, 0, 0, -t111 * t76, t111 * t74, t87 * t118 - t26, 0.2e1 * t91 * t101, t33 - t132, -t29 + t134, -t99, -t6 * qJD(5) - t76 * t139 - t142, t37 * qJD(3) + t4 * qJD(5) - t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, t124, 0, 0, 0, 0, 0, -t114, t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t81 * qJD(5), 0, 0, 0, -pkin(4) * t137, -pkin(4) * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, -t95, t129 + t136, t18, -t125, -pkin(7) * t136 - t102, pkin(7) * t137 - t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t96, t91 * t118 - t128, t74 * t114 + t27, t67, t36 * qJD(3) + t6 * qJD(4) - t140, -t4 * qJD(4) + t74 * t139 - t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t93 * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t95, -t129, t28, t125, t102, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
