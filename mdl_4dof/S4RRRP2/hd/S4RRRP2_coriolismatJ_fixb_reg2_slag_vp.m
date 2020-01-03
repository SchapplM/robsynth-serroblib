% Calculate inertial parameters regressor of coriolis matrix for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:13:11
% EndTime: 2019-12-31 17:13:13
% DurationCPUTime: 0.83s
% Computational Cost: add. (724->147), mult. (1566->170), div. (0->0), fcn. (1096->4), ass. (0->114)
t92 = cos(qJ(2));
t145 = t92 * pkin(1);
t78 = -pkin(2) - t145;
t151 = t78 / 0.2e1 - pkin(2) / 0.2e1;
t118 = qJD(1) + qJD(2);
t89 = sin(qJ(3));
t87 = t89 ^ 2;
t91 = cos(qJ(3));
t88 = t91 ^ 2;
t74 = t88 + t87;
t153 = t118 * t74;
t75 = t88 - t87;
t152 = t118 * t75;
t148 = pkin(3) * t89;
t90 = sin(qJ(2));
t147 = t90 * pkin(1);
t146 = t91 * pkin(3);
t77 = pkin(6) + t147;
t129 = qJ(4) + t77;
t46 = t129 * t89;
t143 = t46 * t89;
t47 = t129 * t91;
t142 = t47 * t91;
t79 = -pkin(2) - t146;
t62 = t79 - t145;
t53 = t62 * t89;
t141 = t62 * t91;
t137 = pkin(6) + qJ(4);
t63 = t137 * t89;
t140 = t63 * t89;
t64 = t137 * t91;
t139 = t64 * t91;
t68 = t79 * t89;
t138 = t79 * t91;
t106 = t74 * t92;
t45 = pkin(1) * t106;
t123 = t45 * qJD(1);
t38 = t45 * qJD(2);
t136 = t123 + t38;
t134 = pkin(1) * qJD(1);
t113 = t90 * t134;
t71 = t89 * t113;
t133 = pkin(1) * qJD(2);
t116 = t90 * t133;
t73 = t89 * t116;
t135 = t71 + t73;
t132 = pkin(2) * qJD(2);
t3 = pkin(3) * t53;
t131 = t3 * qJD(1);
t15 = t142 + t143;
t8 = (t15 * t92 + t62 * t90) * pkin(1);
t130 = t8 * qJD(1);
t128 = qJD(1) * t78;
t127 = t15 * qJD(1);
t16 = (t77 * t106 + t78 * t90) * pkin(1);
t126 = t16 * qJD(1);
t117 = t89 * t146;
t27 = -t53 + t117;
t125 = t27 * qJD(1);
t86 = t87 * pkin(3);
t36 = t86 + t141;
t124 = t36 * qJD(1);
t122 = t47 * qJD(3);
t121 = t64 * qJD(3);
t84 = t89 * qJD(3);
t120 = t89 * qJD(4);
t85 = t91 * qJD(3);
t119 = t91 * qJD(4);
t115 = pkin(3) * t85;
t114 = pkin(3) * t120;
t111 = t147 / 0.2e1;
t110 = -t145 / 0.2e1;
t109 = t89 * t128;
t108 = t91 * t128;
t107 = -t68 / 0.2e1 - t53 / 0.2e1;
t105 = pkin(1) * t118;
t104 = t91 * t116;
t103 = t89 * t110;
t102 = t91 * t110;
t60 = t118 * t89;
t61 = t118 * t91;
t18 = t139 + t140;
t95 = t110 - t79 / 0.2e1 - t62 / 0.2e1;
t1 = t95 * t148;
t9 = pkin(3) * t68;
t101 = -t1 * qJD(1) + t9 * qJD(2);
t4 = t111 + (-t64 / 0.2e1 - t47 / 0.2e1) * t91 + (-t63 / 0.2e1 - t46 / 0.2e1) * t89;
t100 = -t4 * qJD(1) + t18 * qJD(2);
t72 = t91 * t113;
t99 = -t72 - t104;
t10 = (t110 + t146) * t89 + t107;
t37 = -t68 + t117;
t98 = -t10 * qJD(1) - t37 * qJD(2);
t12 = t95 * t91 - t86;
t52 = t86 + t138;
t97 = -t12 * qJD(1) + t52 * qJD(2);
t96 = t110 - t151;
t23 = t96 * t89;
t94 = t23 * qJD(1) + t89 * t132;
t24 = t96 * t91;
t93 = t24 * qJD(1) + t91 * t132;
t83 = pkin(3) * t84;
t76 = t89 * t85;
t66 = t75 * qJD(3);
t65 = t74 * qJD(4);
t54 = pkin(3) * t60;
t42 = t89 * t61;
t26 = t151 * t91 + t102;
t25 = t151 * t89 + t103;
t13 = t86 + t138 / 0.2e1 + t141 / 0.2e1 + t102;
t11 = (t110 - t146) * t89 - t107;
t5 = t139 / 0.2e1 + t142 / 0.2e1 + t140 / 0.2e1 + t143 / 0.2e1 + t111;
t2 = pkin(3) * t103 + (t62 + t79) * t148 / 0.2e1;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t92 * t133, 0, 0, t76, t66, 0, -t76, 0, 0, t78 * t84 - t104, t78 * t85 + t73, t38, t16 * qJD(2), t76, t66, 0, -t76, 0, 0, -t27 * qJD(3) - t104, t36 * qJD(3) + t73, t38 + t65, t8 * qJD(2) + t3 * qJD(3) + t15 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90 * t105, -t92 * t105, 0, 0, t76, t66, 0, -t76, 0, 0, t25 * qJD(3) + t99, t26 * qJD(3) + t135, t136, t126 + (-pkin(2) * t90 + pkin(6) * t106) * t133, t76, t66, 0, -t76, 0, 0, t11 * qJD(3) + t99, t13 * qJD(3) + t135, t136 + t65, t130 + t2 * qJD(3) + t5 * qJD(4) + (t18 * t92 + t79 * t90) * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t152, t85, -t42, -t84, 0, t25 * qJD(2) - t77 * t85 + t109, t26 * qJD(2) + t77 * t84 + t108, 0, 0, t42, t152, t85, -t42, -t84, 0, t11 * qJD(2) - t122 - t125, t13 * qJD(2) + t46 * qJD(3) + t124, -t115, -pkin(3) * t122 + t2 * qJD(2) + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t5 * qJD(2) + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t92 * t134, 0, 0, t76, t66, 0, -t76, 0, 0, -t23 * qJD(3) + t72, -t24 * qJD(3) - t71, -t123, -t126, t76, t66, 0, -t76, 0, 0, -t10 * qJD(3) + t72, -t12 * qJD(3) - t71, -t123 + t65, -t1 * qJD(3) - t4 * qJD(4) - t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t66, 0, -t76, 0, 0, -pkin(2) * t84, -pkin(2) * t85, 0, 0, t76, t66, 0, -t76, 0, 0, -t37 * qJD(3), t52 * qJD(3), t65, t9 * qJD(3) + t18 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t152, t85, -t42, -t84, 0, -pkin(6) * t85 - t94, pkin(6) * t84 - t93, 0, 0, t42, t152, t85, -t42, -t84, 0, t98 - t121, t63 * qJD(3) + t97, -t115, -pkin(3) * t121 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t152, 0, t42, 0, 0, t23 * qJD(2) - t109, t24 * qJD(2) - t108, 0, 0, -t42, -t152, 0, t42, 0, 0, t10 * qJD(2) - t120 + t125, t12 * qJD(2) - t119 - t124, 0, t1 * qJD(2) - t114 - t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t152, 0, t42, 0, 0, t94, t93, 0, 0, -t42, -t152, 0, t42, 0, 0, -t98 - t120, -t97 - t119, 0, -t101 - t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t61, 0, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t85, -t153, t4 * qJD(2) - t127 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t85, -t153, -t100 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t61, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
