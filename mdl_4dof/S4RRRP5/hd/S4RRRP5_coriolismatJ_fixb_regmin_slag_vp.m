% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRP5
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
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP5_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:07
% EndTime: 2019-12-31 17:17:09
% DurationCPUTime: 0.68s
% Computational Cost: add. (1012->117), mult. (2008->149), div. (0->0), fcn. (1907->4), ass. (0->95)
t78 = qJD(2) + qJD(3);
t135 = -pkin(6) - pkin(5);
t82 = cos(qJ(2));
t68 = t135 * t82;
t80 = sin(qJ(3));
t131 = cos(qJ(3));
t81 = sin(qJ(2));
t91 = t131 * t81;
t90 = -t135 * t91 - t80 * t68;
t144 = t78 * t90;
t126 = t80 * t81;
t84 = t135 * t126 - t131 * t68;
t143 = t84 * pkin(3);
t142 = qJ(4) / 0.2e1;
t140 = t84 * qJD(4);
t139 = t78 * t84;
t65 = t80 * t82 + t91;
t62 = t65 ^ 2;
t138 = t84 / 0.2e1;
t133 = t80 * pkin(2);
t72 = qJ(4) + t133;
t137 = -t72 / 0.2e1;
t136 = t80 / 0.2e1;
t63 = -t131 * t82 + t126;
t134 = t63 * pkin(3);
t132 = t81 * pkin(2);
t75 = -t82 * pkin(2) - pkin(1);
t117 = t65 * qJ(4);
t89 = -t117 + t134;
t24 = t75 + t89;
t130 = t24 * t63;
t129 = t24 * t65;
t128 = t72 * t65;
t96 = t131 * pkin(2);
t74 = -t96 - pkin(3);
t127 = t74 * t63;
t121 = pkin(2) * qJD(3);
t32 = t65 * pkin(3) + t63 * qJ(4);
t26 = t32 + t132;
t3 = t24 * t26;
t120 = t3 * qJD(1);
t4 = t24 * t32;
t119 = t4 * qJD(1);
t8 = (t137 + t133 / 0.2e1 + t142) * t65 + (-t96 / 0.2e1 - t74 / 0.2e1 - pkin(3) / 0.2e1) * t63;
t116 = t8 * qJD(1);
t9 = t26 * t63 + t129;
t115 = t9 * qJD(1);
t114 = qJD(1) * t65;
t113 = qJD(1) * t75;
t112 = qJD(1) * t82;
t111 = qJD(3) * t75;
t10 = -t26 * t65 + t130;
t110 = t10 * qJD(1);
t11 = t32 * t63 + t129;
t109 = t11 * qJD(1);
t12 = -t32 * t65 + t130;
t108 = t12 * qJD(1);
t18 = t63 ^ 2 - t62;
t107 = t18 * qJD(1);
t21 = t63 * t132 + t75 * t65;
t106 = t21 * qJD(1);
t22 = t65 * t132 - t75 * t63;
t105 = t22 * qJD(1);
t104 = t62 * qJD(1);
t57 = t63 * qJD(4);
t71 = -t81 ^ 2 + t82 ^ 2;
t103 = t71 * qJD(1);
t102 = t81 * qJD(2);
t101 = t82 * qJD(2);
t77 = qJD(3) * t96;
t100 = t77 + qJD(4);
t99 = pkin(1) * t81 * qJD(1);
t98 = pkin(1) * t112;
t97 = t80 * t121;
t95 = t24 * t114;
t43 = t63 * t114;
t94 = t63 * t113;
t93 = t65 * t113;
t92 = t81 * t112;
t31 = t78 * t65;
t83 = (t131 * t138 + t136 * t90) * pkin(2) + t90 * t137 + t74 * t138;
t85 = t143 / 0.2e1 + t90 * t142;
t2 = t83 + t85;
t50 = (t131 * t72 + t74 * t80) * pkin(2);
t88 = t2 * qJD(1) + t50 * qJD(2);
t79 = qJ(4) * qJD(4);
t76 = qJD(2) * t96;
t70 = t78 * qJ(4);
t69 = t72 * qJD(4);
t30 = t78 * t63;
t23 = qJD(2) * t133;
t17 = t78 * t133;
t7 = -t128 / 0.2e1 - t127 / 0.2e1 - t117 / 0.2e1 + t134 / 0.2e1 + (-t131 * t63 / 0.2e1 + t65 * t136) * pkin(2);
t1 = t83 - t85;
t5 = [0, 0, 0, t81 * t101, t71 * qJD(2), 0, 0, 0, -pkin(1) * t102, -pkin(1) * t101, -t63 * t31, t78 * t18, 0, 0, 0, t21 * qJD(2) + t65 * t111, t22 * qJD(2) - t63 * t111, t9 * qJD(2) + t11 * qJD(3) - t65 * t57, 0, t10 * qJD(2) + t12 * qJD(3) + t62 * qJD(4), t3 * qJD(2) + t4 * qJD(3) - qJD(4) * t129; 0, 0, 0, t92, t103, t101, -t102, 0, -pkin(5) * t101 - t99, pkin(5) * t102 - t98, -t43, t107, -t30, -t31, 0, t106 - t139, t144 + t105, t115 - t139, (-t127 - t128) * qJD(2) + t7 * qJD(3) - t57, -t144 + t110, t120 + (-t72 * t90 + t74 * t84) * qJD(2) + t1 * qJD(3) + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, t107, -t30, -t31, 0, t93 - t139, t144 - t94, t109 - t139, t7 * qJD(2) + t89 * qJD(3) - t57, -t144 + t108, t119 + t1 * qJD(2) + (-qJ(4) * t90 - t143) * qJD(3) + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t30, t104, -t95 + t139; 0, 0, 0, -t92, -t103, 0, 0, 0, t99, t98, t43, -t107, 0, 0, 0, -t106, -t105, -t115, qJD(3) * t8, -t110, qJD(3) * t2 - t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t77, -t97, 0, t100, qJD(3) * t50 + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t76 - t77, -t17, t116, t100 + t76, (-pkin(3) * t80 + t131 * qJ(4)) * t121 + t69 + t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t78 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t107, 0, 0, 0, -t93, t94, -t109, -t8 * qJD(2), -t108, -qJD(2) * t2 - t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t76, t23, -t116, qJD(4) - t76, t79 - t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, -t104, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -qJ(4) * qJD(3) - qJD(2) * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
