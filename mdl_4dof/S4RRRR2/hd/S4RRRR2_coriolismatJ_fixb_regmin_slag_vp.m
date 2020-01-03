% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:20
% EndTime: 2019-12-31 17:23:22
% DurationCPUTime: 0.66s
% Computational Cost: add. (808->117), mult. (1756->158), div. (0->0), fcn. (1664->6), ass. (0->107)
t112 = qJD(3) + qJD(4);
t85 = cos(qJ(2));
t141 = t85 * pkin(1);
t76 = -pkin(2) - t141;
t154 = pkin(2) / 0.2e1 - t76 / 0.2e1;
t113 = qJD(1) + qJD(2);
t139 = cos(qJ(4));
t84 = cos(qJ(3));
t102 = t139 * t84;
t81 = sin(qJ(4));
t82 = sin(qJ(3));
t134 = t81 * t82;
t54 = -t102 + t134;
t103 = t139 * t82;
t133 = t81 * t84;
t56 = -t103 - t133;
t18 = t54 ^ 2 - t56 ^ 2;
t153 = t113 * t18;
t70 = -t82 ^ 2 + t84 ^ 2;
t152 = t113 * t70;
t144 = pkin(6) + pkin(7);
t65 = t144 * t82;
t66 = t144 * t84;
t151 = t112 * (-t139 * t66 + t81 * t65);
t150 = t112 * (t139 * t65 + t81 * t66);
t83 = sin(qJ(2));
t75 = t83 * pkin(1) + pkin(6);
t140 = pkin(7) + t75;
t48 = t140 * t82;
t49 = t140 * t84;
t149 = t112 * (-t139 * t49 + t81 * t48);
t148 = t112 * (t139 * t48 + t81 * t49);
t109 = -t141 / 0.2e1;
t147 = t109 - t154;
t143 = pkin(3) * t81;
t142 = pkin(3) * t82;
t77 = -t84 * pkin(3) - pkin(2);
t63 = t77 - t141;
t138 = t63 * t54;
t137 = t63 * t56;
t136 = t77 * t54;
t135 = t77 * t56;
t132 = pkin(1) * qJD(1);
t131 = pkin(1) * qJD(2);
t130 = pkin(2) * qJD(2);
t129 = qJD(1) * t63;
t128 = qJD(1) * t76;
t127 = qJD(2) * t77;
t126 = qJD(4) * t63;
t125 = qJD(4) * t77;
t45 = t54 * t142;
t16 = t45 - t137;
t120 = t16 * qJD(1);
t46 = t56 * t142;
t17 = -t46 - t138;
t119 = t17 * qJD(1);
t114 = t82 * qJD(3);
t80 = t84 * qJD(3);
t111 = t83 * t131;
t110 = t83 * t132;
t108 = t54 * t129;
t107 = t56 * t129;
t106 = t82 * t128;
t105 = t84 * t128;
t104 = t77 / 0.2e1 + t63 / 0.2e1;
t101 = t139 * qJD(3);
t100 = t139 * qJD(4);
t99 = pkin(1) * t113;
t98 = t54 * t110;
t97 = t56 * t110;
t96 = t82 * t110;
t26 = t112 * t56;
t95 = t83 * t99;
t23 = t45 - t135;
t87 = (-t133 / 0.2e1 - t103 / 0.2e1) * t141;
t8 = t104 * t56 + t87;
t4 = -t45 + t8;
t94 = t4 * qJD(1) - t23 * qJD(2);
t24 = -t46 - t136;
t86 = (-t102 / 0.2e1 + t134 / 0.2e1) * t141;
t9 = t104 * t54 + t86;
t5 = t46 + t9;
t93 = t5 * qJD(1) - t24 * qJD(2);
t92 = t109 + t154;
t31 = t92 * t82;
t91 = t31 * qJD(1) + t82 * t130;
t32 = t92 * t84;
t90 = t32 * qJD(1) + t84 * t130;
t89 = t9 * qJD(1) + t54 * t127;
t88 = t8 * qJD(1) + t56 * t127;
t10 = -t136 / 0.2e1 - t138 / 0.2e1 + t86;
t11 = -t135 / 0.2e1 - t137 / 0.2e1 + t87;
t71 = t82 * t80;
t69 = t82 * t111;
t64 = t70 * qJD(3);
t47 = t113 * t84 * t82;
t44 = t56 * t111;
t43 = t54 * t111;
t34 = t147 * t84;
t33 = t147 * t82;
t25 = t112 * t54;
t13 = t113 * t56 * t54;
t12 = t54 * t26;
t7 = -t46 + t10;
t6 = t45 + t11;
t1 = t112 * t18;
t2 = [0, 0, 0, 0, -t111, -t85 * t131, t71, t64, 0, 0, 0, -t84 * t111 + t76 * t114, t76 * t80 + t69, t12, t1, 0, 0, 0, t16 * qJD(3) - t56 * t126 + t43, t17 * qJD(3) - t54 * t126 - t44; 0, 0, 0, 0, -t95, -t85 * t99, t71, t64, 0, 0, 0, t33 * qJD(3) - t84 * t95, t34 * qJD(3) + t69 + t96, t12, t1, 0, 0, 0, t6 * qJD(3) + t11 * qJD(4) + t43 + t98, t7 * qJD(3) + t10 * qJD(4) - t44 - t97; 0, 0, 0, 0, 0, 0, t47, t152, t80, -t114, 0, t33 * qJD(2) - t75 * t80 + t106, t34 * qJD(2) + t75 * t114 + t105, t13, t153, -t25, t26, 0, t6 * qJD(2) + t120 + t149, t7 * qJD(2) + t119 + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t153, -t25, t26, 0, t11 * qJD(2) - t107 + t149, t10 * qJD(2) - t108 + t148; 0, 0, 0, 0, t110, t85 * t132, t71, t64, 0, 0, 0, -t31 * qJD(3) + t84 * t110, -t32 * qJD(3) - t96, t12, t1, 0, 0, 0, -t4 * qJD(3) - t8 * qJD(4) - t98, -t5 * qJD(3) - t9 * qJD(4) + t97; 0, 0, 0, 0, 0, 0, t71, t64, 0, 0, 0, -pkin(2) * t114, -pkin(2) * t80, t12, t1, 0, 0, 0, t23 * qJD(3) - t56 * t125, t24 * qJD(3) - t54 * t125; 0, 0, 0, 0, 0, 0, t47, t152, t80, -t114, 0, -pkin(6) * t80 - t91, pkin(6) * t114 - t90, t13, t153, -t25, t26, 0, -t94 + t151, -t93 + t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, t153, -t25, t26, 0, -t88 + t151, -t89 + t150; 0, 0, 0, 0, 0, 0, -t47, -t152, 0, 0, 0, t31 * qJD(2) - t106, t32 * qJD(2) - t105, -t13, -t153, 0, 0, 0, t4 * qJD(2) - t120, t5 * qJD(2) - t119; 0, 0, 0, 0, 0, 0, -t47, -t152, 0, 0, 0, t91, t90, -t13, -t153, 0, 0, 0, t94, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t143, -pkin(3) * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112 * t143, (-t101 - t100) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t153, 0, 0, 0, t8 * qJD(2) + t107, t9 * qJD(2) + t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t153, 0, 0, 0, t88, t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t143, pkin(3) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t2;
