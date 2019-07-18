% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:49
% EndTime: 2019-07-18 13:28:54
% DurationCPUTime: 1.27s
% Computational Cost: add. (949->164), mult. (2726->270), div. (0->0), fcn. (2200->8), ass. (0->103)
t65 = sin(qJ(3));
t112 = qJD(2) * t65;
t68 = cos(qJ(3));
t131 = cos(qJ(4));
t97 = qJD(2) * t131;
t57 = t68 * t97;
t64 = sin(qJ(4));
t137 = -t64 * t112 + t57;
t60 = qJD(3) + qJD(4);
t18 = t137 * t60;
t38 = qJD(5) - t137;
t136 = qJD(5) - t38;
t45 = t131 * t65 + t64 * t68;
t135 = t60 * t45;
t19 = t135 * qJD(2);
t108 = qJD(3) * t65;
t66 = sin(qJ(2));
t114 = qJD(1) * t66;
t130 = t19 * t68;
t134 = (t38 * t108 - t130) * pkin(2) - t38 * t114;
t78 = t131 * t68 - t64 * t65;
t20 = t60 * t78;
t100 = t65 * t114;
t49 = qJD(3) * pkin(2) - t100;
t113 = qJD(1) * t68;
t103 = t64 * t113;
t92 = t66 * t103;
t27 = -t131 * t49 + t92;
t69 = cos(qJ(2));
t76 = t78 * t69;
t33 = qJD(1) * t76;
t105 = t69 * qJD(1);
t110 = qJD(2) * t68;
t50 = -pkin(2) * t110 - t105;
t109 = qJD(2) * t69;
t119 = t66 * t68;
t34 = (-qJD(3) * t119 - t65 * t109) * qJD(1);
t89 = t131 * t114;
t96 = t131 * qJD(4);
t6 = -qJD(4) * t92 + t57 * t105 - t89 * t108 + t64 * t34 + t49 * t96;
t93 = t64 * t100;
t117 = -qJD(3) * t93 - t131 * t34;
t121 = t64 * t49;
t7 = (t64 * t109 + t66 * t96) * t113 + qJD(4) * t121 + t117;
t133 = (qJD(5) * t50 + t6) * t78 + t27 * t20 + t7 * t45 + (qJD(5) * t68 * pkin(2) + t33) * t38;
t63 = sin(qJ(5));
t41 = -t64 * t110 - t65 * t97;
t67 = cos(qJ(5));
t84 = t67 * t41 - t63 * t60;
t9 = -t84 * qJD(5) + t63 * t18;
t106 = qJD(5) * t67;
t107 = qJD(5) * t63;
t8 = t60 * t106 + t41 * t107 + t67 * t18;
t132 = t8 * t63;
t24 = -t63 * t41 - t67 * t60;
t129 = t24 * t38;
t128 = t84 * t38;
t127 = t27 * t45;
t126 = t38 * t41;
t125 = t41 * t137;
t124 = t50 * t41;
t122 = t63 * t19;
t120 = t65 * t68;
t118 = t67 * t19;
t116 = t65 ^ 2 - t68 ^ 2;
t70 = qJD(3) ^ 2;
t71 = qJD(2) ^ 2;
t115 = t70 + t71;
t111 = qJD(2) * t66;
t104 = qJD(2) * qJD(3);
t99 = 0.2e1 * t104;
t98 = t27 * (-t38 - t137);
t94 = t38 * t67;
t91 = -0.2e1 * t69 * t104;
t51 = t68 * t89;
t28 = t51 + t121;
t16 = t67 * t28 + t63 * t50;
t90 = t27 * t106 - t16 * t41 + t7 * t63;
t32 = t45 * t114;
t88 = -t137 * t27 - t32 * t38;
t87 = t63 * t28 - t67 * t50;
t36 = t78 * t66;
t86 = t67 * t36 - t69 * t63;
t85 = -t63 * t36 - t69 * t67;
t83 = qJD(5) * t64 + t112;
t81 = t27 * t107 - t41 * t87 - t7 * t67;
t79 = -t45 * t107 + t67 * t20;
t77 = t45 * t69;
t72 = -t137 * t50 - t6;
t43 = (pkin(2) * t108 + t114) * qJD(2);
t37 = t67 * t43;
t35 = t45 * t66;
t31 = qJD(1) * t77;
t30 = -t51 + t93;
t14 = -t137 ^ 2 + t41 ^ 2;
t13 = -t41 * t60 - t19;
t11 = qJD(2) * t77 + t20 * t66;
t10 = qJD(2) * t76 - t135 * t66;
t4 = t38 * t94 - t41 * t84 + t122;
t3 = -t38 ^ 2 * t63 - t24 * t41 + t118;
t2 = -t84 * t94 + t132;
t1 = (t8 - t129) * t67 + (-t9 + t128) * t63;
t5 = [0, 0, -t71 * t66, -t71 * t69, 0, 0, 0, 0, 0, -t115 * t119 + t65 * t91, t115 * t66 * t65 + t68 * t91, 0, 0, 0, 0, 0, -t11 * t60 - t111 * t137 - t69 * t19, -t10 * t60 - t41 * t111 - t69 * t18, 0, 0, 0, 0, 0, (-t86 * qJD(5) - t63 * t10 + t67 * t111) * t38 + t85 * t19 + t11 * t24 + t35 * t9, -(t85 * qJD(5) + t67 * t10 + t63 * t111) * t38 - t86 * t19 - t11 * t84 + t35 * t8; 0, 0, 0, 0, t99 * t120, -t116 * t99, t70 * t68, -t70 * t65, 0, 0, 0, t18 * t45 - t41 * t20, t135 * t41 + t137 * t20 + t18 * t78 - t45 * t19, t20 * t60, -t135 * t60, 0, t137 * t114 + t50 * t135 + t31 * t60 - t43 * t78 + (-t108 * t137 - t130) * pkin(2), t41 * t114 + t50 * t20 + t33 * t60 + t43 * t45 + (-t41 * t108 - t18 * t68) * pkin(2), t8 * t67 * t45 - t79 * t84, (-t24 * t67 + t63 * t84) * t20 + (-t132 - t67 * t9 + (t24 * t63 + t67 * t84) * qJD(5)) * t45, t45 * t118 - t135 * t84 + t79 * t38 - t78 * t8, -t45 * t122 - t24 * t135 + t9 * t78 + (-t45 * t106 - t63 * t20) * t38, t135 * t38 - t19 * t78, -t87 * t135 - t31 * t24 - t37 * t78 + ((t28 * t78 + t127) * qJD(5) + t134) * t67 + t133 * t63, -t16 * t135 + t31 * t84 + t133 * t67 + ((-qJD(5) * t28 + t43) * t78 - qJD(5) * t127 - t134) * t63; 0, 0, 0, 0, -t71 * t120, t116 * t71, 0, 0, 0, 0, 0, t125, t14, 0, t13, 0, -t30 * t60 + t124 + (pkin(2) * t137 * t65 - t69 * t103) * qJD(2) + (-t51 + (-pkin(2) * t60 - t49) * t64) * qJD(4) - t117, -t32 * t60 + (t41 * t112 - t60 * t96) * pkin(2) + t72, t2, t1, t4, t3, t126, t30 * t24 + t88 * t63 + (-t131 * t9 + (qJD(4) * t24 - t122) * t64 + (-t63 * t96 - t83 * t67) * t38) * pkin(2) + t81, -t30 * t84 + t88 * t67 + (-t131 * t8 + (-qJD(4) * t84 - t118) * t64 + (t83 * t63 - t67 * t96) * t38) * pkin(2) + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t14, 0, t13, 0, t28 * t60 + t124 - t7, -t27 * t60 + t72, t2, t1, t4, t3, t126, -t28 * t24 + t63 * t98 + t81, t28 * t84 + t67 * t98 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84 * t24, -t24 ^ 2 + t84 ^ 2, t8 + t129, -t9 - t128, t19, -t136 * t16 + t27 * t84 - t63 * t6 + t37, t136 * t87 + t27 * t24 - t63 * t43 - t67 * t6;];
tauc_reg  = t5;
