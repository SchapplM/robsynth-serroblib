% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR8_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:08:26
% EndTime: 2019-12-31 17:08:28
% DurationCPUTime: 0.93s
% Computational Cost: add. (1070->107), mult. (2017->141), div. (0->0), fcn. (1858->4), ass. (0->93)
t108 = qJD(2) - qJD(4);
t152 = pkin(5) - pkin(6);
t87 = cos(qJ(2));
t75 = t152 * t87;
t84 = sin(qJ(4));
t86 = cos(qJ(4));
t85 = sin(qJ(2));
t96 = t152 * t85;
t150 = t86 * t75 + t84 * t96;
t154 = t108 * t150;
t146 = t86 * t96;
t142 = t146 / 0.2e1;
t64 = -t87 * t84 + t85 * t86;
t144 = t108 * t64;
t62 = t85 * t84 + t87 * t86;
t151 = t62 * t144;
t18 = t62 ^ 2 - t64 ^ 2;
t149 = t18 * qJD(1);
t137 = t84 * t75;
t148 = -t146 + t137;
t145 = t108 * t62;
t143 = t150 / 0.2e1;
t140 = pkin(2) + pkin(3);
t139 = t150 * t86;
t134 = t85 * qJ(3);
t55 = t140 * t87 + pkin(1) + t134;
t138 = t55 * t62;
t133 = t87 * qJ(3);
t60 = -t140 * t85 + t133;
t5 = t55 * t60;
t135 = t5 * qJD(1);
t132 = qJD(1) * t64;
t131 = qJD(1) * t85;
t130 = qJD(1) * t87;
t129 = qJD(3) * t85;
t128 = qJD(4) * t55;
t33 = t139 / 0.2e1;
t12 = t33 - t139 / 0.2e1;
t127 = t12 * qJD(1);
t126 = t12 * qJD(3);
t13 = -t55 * t64 + t60 * t62;
t125 = t13 * qJD(1);
t14 = t60 * t64 + t138;
t124 = t14 * qJD(1);
t90 = -t87 * pkin(2) - t134;
t69 = -pkin(1) + t90;
t72 = t85 * pkin(2) - t133;
t37 = t69 * t87 + t72 * t85;
t121 = t37 * qJD(1);
t38 = -t69 * t85 + t72 * t87;
t120 = t38 * qJD(1);
t40 = t142 - t146 / 0.2e1;
t119 = t40 * qJD(1);
t83 = t85 ^ 2;
t76 = t87 ^ 2 - t83;
t118 = t76 * qJD(1);
t117 = t76 * qJD(2);
t116 = t83 * qJD(1);
t115 = t84 * qJD(2);
t114 = t84 * qJD(3);
t113 = t85 * qJD(2);
t112 = t86 * qJD(2);
t111 = t86 * qJD(3);
t81 = t87 * qJD(2);
t110 = t87 * qJD(3);
t109 = qJD(2) * qJ(3);
t107 = pkin(1) * t131;
t106 = pkin(1) * t130;
t105 = pkin(5) * t113;
t104 = pkin(5) * t81;
t103 = qJD(1) * t138;
t102 = t55 * t132;
t101 = t55 * t131;
t36 = t62 * t132;
t100 = t69 * t72 * qJD(1);
t99 = t69 * t131;
t98 = t62 * t131;
t97 = t64 * t131;
t67 = t84 * qJ(3) + t140 * t86;
t68 = t86 * qJ(3) - t140 * t84;
t35 = t67 * t84 + t68 * t86;
t6 = (t143 - t150 / 0.2e1) * t86;
t89 = -t6 * qJD(1) + t35 * qJD(2);
t88 = qJD(2) * t90 + t110;
t78 = t85 * t81;
t77 = t85 * t130;
t71 = t108 * t86;
t70 = t108 * t84;
t19 = -t137 + 0.2e1 * t142;
t15 = -t86 * t62 + t84 * t64;
t11 = t12 * qJD(4);
t7 = t86 * t143 + t148 * t84 + t33;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t117, 0, -t78, 0, 0, -pkin(1) * t113, -pkin(1) * t81, 0, 0, t78, 0, -t117, 0, 0, -t78, -t38 * qJD(2) + t110 * t85, 0, -t37 * qJD(2) + t83 * qJD(3), (qJD(2) * t72 - t129) * t69, t151, -t108 * t18, 0, -t151, 0, 0, t13 * qJD(2) + t128 * t64 + t129 * t62, t14 * qJD(2) - t128 * t62 + t129 * t64, 0, t5 * qJD(2) + t129 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t118, t81, -t77, -t113, 0, -t104 - t107, t105 - t106, 0, 0, t77, t81, -t118, 0, t113, -t77, -t104 - t120, t88, -t105 - t121, pkin(5) * t88 + t100, t36, -t149, -t145, -t36, -t144, 0, t125 - t154, qJD(2) * t148 + t19 * qJD(4) + t124, (t67 * t62 + t68 * t64) * qJD(2) + t15 * qJD(3), t135 + (t148 * t68 - t150 * t67) * qJD(2) + t7 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t81, t116, -t99 + t104, 0, 0, 0, 0, 0, 0, t98, t97, t15 * qJD(2), t7 * qJD(2) + t101 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t149, t145, t36, t144, 0, t102 + t154, t19 * qJD(2) + qJD(4) * t148 - t103, 0, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t118, 0, t77, 0, 0, t107, t106, 0, 0, -t77, 0, t118, 0, 0, t77, t120, 0, t121, -t100, -t36, t149, 0, t36, 0, 0, -t125, t40 * qJD(4) - t124, 0, -t6 * qJD(3) - t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), 0, 0, 0, 0, 0, 0, t68 * qJD(4) + t114, -t67 * qJD(4) + t111, 0, t35 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t109, 0, 0, 0, 0, 0, 0, t115, t112, 0, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108 * t68, -t108 * t67 + t119, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, 0, -t116, t99, 0, 0, 0, 0, 0, 0, -t98, -t97, 0, t6 * qJD(2) - t101 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t109, 0, 0, 0, 0, 0, 0, -t70, -t71, 0, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t71, 0, t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t149, 0, -t36, 0, 0, -t102, -t40 * qJD(2) + t103, 0, -t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68 * qJD(2) - t114, t67 * qJD(2) - t111 - t119, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -t112, 0, -t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
