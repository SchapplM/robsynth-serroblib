% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRR4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:08:01
% EndTime: 2019-12-05 17:08:04
% DurationCPUTime: 0.75s
% Computational Cost: add. (822->129), mult. (1782->158), div. (0->0), fcn. (1698->6), ass. (0->107)
t109 = qJD(4) + qJD(5);
t85 = cos(qJ(3));
t141 = t85 * pkin(2);
t75 = -pkin(3) - t141;
t154 = pkin(3) / 0.2e1 - t75 / 0.2e1;
t83 = cos(qJ(5));
t84 = cos(qJ(4));
t132 = t83 * t84;
t80 = sin(qJ(5));
t81 = sin(qJ(4));
t135 = t80 * t81;
t55 = -t132 + t135;
t153 = t109 * t55;
t110 = qJD(2) + qJD(3);
t133 = t83 * t81;
t134 = t80 * t84;
t57 = t133 + t134;
t18 = t55 ^ 2 - t57 ^ 2;
t152 = t110 * t18;
t71 = -t81 ^ 2 + t84 ^ 2;
t151 = t110 * t71;
t143 = pkin(7) + pkin(8);
t66 = t143 * t81;
t67 = t143 * t84;
t150 = t109 * (t80 * t66 - t83 * t67);
t149 = t109 * (t83 * t66 + t80 * t67);
t82 = sin(qJ(3));
t74 = t82 * pkin(2) + pkin(7);
t140 = pkin(8) + t74;
t49 = t140 * t81;
t50 = t140 * t84;
t148 = t109 * (t80 * t49 - t83 * t50);
t147 = t109 * (t83 * t49 + t80 * t50);
t106 = -t141 / 0.2e1;
t146 = t106 - t154;
t142 = pkin(4) * t81;
t76 = -t84 * pkin(4) - pkin(3);
t64 = t76 - t141;
t139 = t64 * t55;
t138 = t64 * t57;
t137 = t76 * t55;
t136 = t76 * t57;
t131 = pkin(2) * qJD(2);
t130 = pkin(2) * qJD(3);
t129 = pkin(3) * qJD(3);
t128 = pkin(4) * qJD(5);
t127 = qJD(4) * pkin(4);
t126 = qJD(2) * t64;
t125 = qJD(2) * t75;
t124 = qJD(3) * t76;
t46 = t55 * t142;
t16 = t46 + t138;
t119 = t16 * qJD(2);
t47 = t57 * t142;
t17 = t47 - t139;
t118 = t17 * qJD(2);
t113 = t55 * qJD(5);
t112 = t57 * qJD(5);
t111 = t81 * qJD(4);
t79 = t84 * qJD(4);
t108 = t82 * t130;
t107 = t82 * t131;
t105 = t55 * t126;
t104 = t57 * t126;
t103 = t81 * t125;
t102 = t84 * t125;
t101 = t76 / 0.2e1 + t64 / 0.2e1;
t100 = pkin(2) * t110;
t99 = pkin(4) * t109;
t98 = t55 * t107;
t97 = t57 * t107;
t96 = t81 * t107;
t27 = t109 * t57;
t95 = t82 * t100;
t23 = t46 + t136;
t87 = (-t134 / 0.2e1 - t133 / 0.2e1) * t141;
t8 = -t101 * t57 + t87;
t4 = -t46 + t8;
t94 = t4 * qJD(2) - t23 * qJD(3);
t24 = t47 - t137;
t86 = (-t132 / 0.2e1 + t135 / 0.2e1) * t141;
t9 = t101 * t55 + t86;
t5 = -t47 + t9;
t93 = t5 * qJD(2) - t24 * qJD(3);
t92 = t106 + t154;
t32 = t92 * t81;
t91 = t32 * qJD(2) + t81 * t129;
t33 = t92 * t84;
t90 = t33 * qJD(2) + t84 * t129;
t89 = t9 * qJD(2) + t55 * t124;
t88 = t8 * qJD(2) - t57 * t124;
t10 = -t137 / 0.2e1 - t139 / 0.2e1 + t86;
t11 = t136 / 0.2e1 + t138 / 0.2e1 + t87;
t72 = t81 * t79;
t70 = t81 * t108;
t65 = t71 * qJD(4);
t48 = t110 * t84 * t81;
t45 = t57 * t108;
t44 = t55 * t108;
t35 = t146 * t84;
t34 = t146 * t81;
t13 = t110 * t57 * t55;
t12 = t55 * t27;
t7 = t47 + t10;
t6 = t46 + t11;
t1 = t109 * t18;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t79, 0, 0, 0, 0, 0, -t27, t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t108, -t85 * t130, t72, t65, 0, 0, 0, -t84 * t108 + t75 * t111, t75 * t79 + t70, -t12, t1, 0, 0, 0, t16 * qJD(4) + t64 * t112 + t44, t17 * qJD(4) - t64 * t113 + t45; 0, 0, 0, 0, 0, -t95, -t85 * t100, t72, t65, 0, 0, 0, t34 * qJD(4) - t84 * t95, t35 * qJD(4) + t70 + t96, -t12, t1, 0, 0, 0, t6 * qJD(4) + t11 * qJD(5) + t44 + t98, t7 * qJD(4) + t10 * qJD(5) + t45 + t97; 0, 0, 0, 0, 0, 0, 0, t48, t151, t79, -t111, 0, t34 * qJD(3) - t74 * t79 + t103, t35 * qJD(3) + t74 * t111 + t102, -t13, t152, -t153, -t27, 0, t6 * qJD(3) + t119 + t148, t7 * qJD(3) + t118 + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t152, -t153, -t27, 0, t11 * qJD(3) + t104 + t148, t10 * qJD(3) - t105 + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t107, t85 * t131, t72, t65, 0, 0, 0, -t32 * qJD(4) + t84 * t107, -t33 * qJD(4) - t96, -t12, t1, 0, 0, 0, -t4 * qJD(4) - t8 * qJD(5) - t98, -t5 * qJD(4) - t9 * qJD(5) - t97; 0, 0, 0, 0, 0, 0, 0, t72, t65, 0, 0, 0, -pkin(3) * t111, -pkin(3) * t79, -t12, t1, 0, 0, 0, t23 * qJD(4) + t76 * t112, t24 * qJD(4) - t76 * t113; 0, 0, 0, 0, 0, 0, 0, t48, t151, t79, -t111, 0, -pkin(7) * t79 - t91, pkin(7) * t111 - t90, -t13, t152, -t153, -t27, 0, -t94 + t150, -t93 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t152, -t153, -t27, 0, -t88 + t150, -t89 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t48, -t151, 0, 0, 0, t32 * qJD(3) - t103, t33 * qJD(3) - t102, t13, -t152, 0, 0, 0, t4 * qJD(3) - t119, t5 * qJD(3) - t118; 0, 0, 0, 0, 0, 0, 0, -t48, -t151, 0, 0, 0, t91, t90, t13, -t152, 0, 0, 0, t94, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 * t128, -t83 * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 * t99, -t83 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t152, 0, 0, 0, t8 * qJD(3) - t104, t9 * qJD(3) + t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t152, 0, 0, 0, t88, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t127, t83 * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
