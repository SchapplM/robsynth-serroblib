% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x22]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPPR6_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:53
% EndTime: 2019-12-31 17:47:56
% DurationCPUTime: 0.93s
% Computational Cost: add. (495->128), mult. (1312->213), div. (0->0), fcn. (1260->6), ass. (0->123)
t152 = pkin(3) + qJ(2);
t75 = sin(pkin(8));
t71 = t75 ^ 2;
t77 = cos(pkin(8));
t73 = t77 ^ 2;
t138 = t71 + t73;
t78 = cos(pkin(7));
t80 = cos(qJ(5));
t140 = t80 * t78;
t76 = sin(pkin(7));
t79 = sin(qJ(5));
t144 = t76 * t79;
t38 = t75 * t140 - t144;
t151 = t38 / 0.2e1;
t141 = t79 * t78;
t143 = t76 * t80;
t37 = t75 * t141 + t143;
t149 = t37 * t38;
t74 = t78 ^ 2;
t148 = t73 * t74;
t147 = t75 * t37;
t146 = t75 * t38;
t145 = t75 * t78;
t142 = t77 * t78;
t91 = -t76 * qJ(3) - pkin(1);
t39 = (-pkin(2) - qJ(4)) * t78 + t91;
t58 = t152 * t76;
t25 = t77 * t39 + t75 * t58;
t139 = t152 * t78;
t72 = t76 ^ 2;
t62 = t72 + t74;
t22 = t76 * pkin(6) + t25;
t28 = (pkin(4) * t77 + pkin(6) * t75) * t78 + t139;
t6 = t79 * t22 - t80 * t28;
t86 = -t75 * t39 + t77 * t58;
t81 = t76 * pkin(4) + t86;
t1 = -t6 * t142 + t37 * t81;
t137 = t1 * qJD(1);
t7 = -t80 * t22 - t79 * t28;
t2 = t7 * t142 + t38 * t81;
t136 = t2 * qJD(1);
t3 = t139 * t78 + (t25 * t75 + t86 * t77) * t76;
t135 = t3 * qJD(1);
t83 = t75 * t86;
t4 = (t25 * t77 - t83) * t76;
t134 = t4 * qJD(1);
t5 = t25 * t142 - t78 * t83;
t133 = t5 * qJD(1);
t9 = t37 ^ 2 - t38 ^ 2;
t132 = t9 * qJD(1);
t131 = qJD(4) * t78;
t130 = qJD(5) * t79;
t129 = qJD(5) * t80;
t106 = t73 / 0.2e1 + 0.1e1 / 0.2e1;
t88 = t78 * t106;
t96 = -t144 / 0.2e1;
t10 = t80 * t88 + (t96 + t151) * t75;
t128 = t10 * qJD(1);
t94 = t143 / 0.2e1;
t11 = t79 * t88 + (t37 / 0.2e1 + t94) * t75;
t127 = t11 * qJD(1);
t14 = ((-t75 * t144 + t140) * t78 + t76 * t37) * t77;
t126 = t14 * qJD(1);
t15 = ((t75 * t143 + t141) * t78 - t76 * t38) * t77;
t125 = t15 * qJD(1);
t93 = t141 / 0.2e1;
t85 = t75 * t93 - t37 / 0.2e1;
t95 = -t143 / 0.2e1;
t17 = (t95 + t85) * t77;
t124 = t17 * qJD(1);
t92 = -t140 / 0.2e1;
t84 = t75 * t92 + t151;
t19 = (t96 + t84) * t77;
t123 = t19 * qJD(1);
t23 = (t73 * t141 + t147) * t76;
t122 = t23 * qJD(1);
t24 = (t73 * t140 + t146) * t76;
t121 = t24 * qJD(1);
t26 = t37 * t145 + t79 * t148;
t120 = t26 * qJD(1);
t27 = t38 * t145 + t80 * t148;
t119 = t27 * qJD(1);
t40 = t138 * t76;
t33 = t78 * t40;
t118 = t33 * qJD(1);
t35 = (t71 / 0.2e1 + t106) * t78;
t117 = t35 * qJD(1);
t116 = t40 * qJD(1);
t41 = t62 * t75;
t115 = t41 * qJD(1);
t42 = t138 * t74;
t114 = t42 * qJD(1);
t43 = t62 * t77;
t113 = t43 * qJD(1);
t57 = t62 * qJ(2);
t112 = t57 * qJD(1);
t111 = t62 * qJD(1);
t110 = t72 * qJD(1);
t109 = t72 * qJD(3);
t108 = t76 * qJD(1);
t107 = t76 * qJD(3);
t105 = qJD(1) * t142;
t104 = t77 * t131;
t103 = qJD(5) * t142;
t102 = t77 * t130;
t101 = t77 * t129;
t100 = qJD(1) * t149;
t99 = t75 * t110;
t98 = t78 * t108;
t97 = t77 * t110;
t90 = t75 * t98;
t89 = t77 * t98;
t87 = qJD(5) + t105;
t82 = t76 * t131 + t109;
t59 = t62 * qJD(2);
t50 = -t78 * pkin(2) + t91;
t47 = t57 * qJD(2);
t36 = (0.1e1 / 0.2e1 - t138 / 0.2e1) * t78;
t18 = (t144 / 0.2e1 + t84) * t77;
t16 = (t85 + t94) * t77;
t13 = t73 * t92 - t146 / 0.2e1 + t75 * t96 + t140 / 0.2e1;
t12 = t73 * t93 + t147 / 0.2e1 + t75 * t95 - t141 / 0.2e1;
t8 = [0, 0, 0, 0, 0, t59, t47, t59, -t78 * t107, t109, -t50 * t107 + t47, t43 * qJD(2) + t82 * t75, -t41 * qJD(2) + t82 * t77, t33 * qJD(3) + t42 * qJD(4), t3 * qJD(2) - t4 * qJD(3) - t5 * qJD(4), -qJD(5) * t149, t9 * qJD(5), t37 * t103, t38 * t103, 0, t14 * qJD(2) + t23 * qJD(3) + t26 * qJD(4) + t2 * qJD(5), -t15 * qJD(2) + t24 * qJD(3) + t27 * qJD(4) - t1 * qJD(5); 0, 0, 0, 0, 0, t111, t112, t111, 0, 0, t112, t113, -t115, 0, t36 * qJD(4) + t135, 0, 0, 0, 0, 0, t13 * qJD(5) + t126, t12 * qJD(5) - t125; 0, 0, 0, 0, 0, 0, 0, 0, -t98, t110, -t50 * t108, t99, t97, t118, -t134, 0, 0, 0, 0, 0, t18 * qJD(5) + t122, t16 * qJD(5) + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t89, t114, t36 * qJD(2) - t133, 0, 0, 0, 0, 0, t120, t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t132, t87 * t37, t87 * t38, 0, t13 * qJD(2) + t18 * qJD(3) + t7 * qJD(5) + t136, t12 * qJD(2) + t16 * qJD(3) + t6 * qJD(5) - t137; 0, 0, 0, 0, 0, -t111, -t112, -t111, 0, 0, -t112 - t107, -t113, t115, 0, -t40 * qJD(3) - t35 * qJD(4) - t135, 0, 0, 0, 0, 0, -t10 * qJD(5) - t126, t11 * qJD(5) + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, 0, 0, 0, -t116, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 - t128, t102 + t127; 0, 0, 0, 0, 0, 0, 0, 0, t98, -t110, (qJD(1) * t50 + qJD(2)) * t76, -t99, -t97, -t118, t40 * qJD(2) + t134, 0, 0, 0, 0, 0, t19 * qJD(5) - t122, t17 * qJD(5) - t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, 0, 0, 0, t116, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129 * t75 + t123, t130 * t75 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, -t89, -t114, t35 * qJD(2) + t133, 0, 0, 0, 0, 0, -t78 * t102 - t120, -t78 * t101 - t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87 * t79, -t87 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, -t132, -t37 * t105, -t38 * t105, 0, t10 * qJD(2) - t19 * qJD(3) + t104 * t79 - t136, -t11 * qJD(2) - t17 * qJD(3) + t104 * t80 + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, -t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, -t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79 * t105, t80 * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t8;
