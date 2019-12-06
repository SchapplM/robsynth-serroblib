% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:59
% EndTime: 2019-12-05 17:40:03
% DurationCPUTime: 0.87s
% Computational Cost: add. (991->94), mult. (1807->122), div. (0->0), fcn. (2082->6), ass. (0->83)
t151 = qJD(4) + qJD(5);
t130 = sin(qJ(4));
t73 = sin(pkin(8));
t74 = cos(pkin(8));
t78 = cos(qJ(4));
t56 = t130 * t74 + t78 * t73;
t76 = sin(qJ(5));
t126 = t76 * t56;
t58 = -t130 * t73 + t78 * t74;
t77 = cos(qJ(5));
t47 = t77 * t58;
t33 = t47 - t126;
t125 = t76 * t58;
t46 = t77 * t56;
t34 = t46 + t125;
t145 = -t33 ^ 2 + t34 ^ 2;
t150 = t145 * qJD(1);
t149 = t33 * qJD(1);
t148 = t33 * qJD(4);
t104 = t34 * qJD(4);
t106 = t34 * qJD(5);
t9 = -t104 - t106;
t75 = -pkin(1) - qJ(3);
t131 = -pkin(6) + t75;
t59 = t131 * t73;
t60 = t131 * t74;
t80 = -t130 * t60 - t78 * t59;
t26 = -t56 * pkin(7) - t80;
t83 = t130 * t59 - t78 * t60;
t79 = -t58 * pkin(7) - t83;
t147 = t151 * (-t77 * t26 - t76 * t79);
t146 = t151 * (t76 * t26 - t77 * t79);
t141 = t34 * t149;
t138 = qJD(3) * t34;
t137 = t33 * qJD(5);
t136 = t34 * qJD(1);
t135 = -t137 - t148;
t89 = -t46 / 0.2e1;
t88 = t47 / 0.2e1;
t132 = pkin(4) * t58;
t69 = t73 ^ 2;
t70 = t74 ^ 2;
t63 = t69 + t70;
t121 = pkin(4) * qJD(5);
t120 = qJD(4) * pkin(4);
t64 = t73 * pkin(3) + qJ(2);
t38 = t56 * pkin(4) + t64;
t6 = -t34 * t132 - t38 * t33;
t117 = t6 * qJD(1);
t7 = -t132 * t33 + t34 * t38;
t116 = t7 * qJD(1);
t15 = 0.2e1 * t88 - t126;
t112 = t15 * qJD(1);
t27 = t56 ^ 2 - t58 ^ 2;
t111 = t27 * qJD(1);
t32 = t88 - t47 / 0.2e1;
t110 = t32 * qJD(1);
t109 = t32 * qJD(5);
t55 = t63 * t75;
t101 = t55 * qJD(1);
t100 = t56 * qJD(1);
t49 = t56 * qJD(4);
t99 = t58 * qJD(1);
t52 = t58 * qJD(4);
t86 = -t69 / 0.2e1 - t70 / 0.2e1;
t62 = -0.1e1 / 0.2e1 + t86;
t98 = t62 * qJD(1);
t97 = t63 * qJD(1);
t96 = t73 * qJD(1);
t95 = t74 * qJD(1);
t92 = t38 * t136;
t91 = t38 * t149;
t90 = t56 * t99;
t85 = pkin(4) * t151;
t84 = qJD(1) * t64 + qJD(3);
t31 = t89 + t46 / 0.2e1;
t82 = t31 * qJD(2);
t81 = t15 * qJD(5) + t148;
t72 = qJ(2) * qJD(2);
t71 = qJD(1) * qJ(2);
t61 = 0.1e1 / 0.2e1 + t86;
t14 = 0.2e1 * t89 - t125;
t1 = [0, 0, 0, 0, qJD(2), t72, qJD(2) * t73, qJD(2) * t74, t63 * qJD(3), -t55 * qJD(3) + t72, -t56 * t52, t27 * qJD(4), 0, 0, 0, qJD(2) * t56 + t64 * t52, qJD(2) * t58 - t64 * t49, t9 * t33, t151 * t145, 0, 0, 0, t34 * qJD(2) - t6 * qJD(4) + t137 * t38, qJD(2) * t33 - t7 * qJD(4) - t38 * t106; 0, 0, 0, 0, qJD(1), t71, t96, t95, 0, t61 * qJD(3) + t71, 0, 0, 0, 0, 0, t100, t99, 0, 0, 0, 0, 0, t136, t149; 0, 0, 0, 0, 0, 0, 0, 0, t97, t61 * qJD(2) - t101, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90, t111, -t49, -t52, 0, t80 * qJD(4) + t64 * t99, t83 * qJD(4) - t64 * t100, -t141, t150, t9, -t81, 0, -t117 + t147, -t116 + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141, t150, t9, -t15 * qJD(4) - t137, 0, t32 * qJD(3) + t147 + t91, -t92 + t146; 0, 0, 0, 0, -qJD(1), -t71, -t96, -t95, 0, t62 * qJD(3) - t71, 0, 0, 0, 0, 0, -t100, -t99, 0, 0, 0, 0, 0, -t136, -t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t52, 0, 0, 0, 0, 0, t14 * qJD(5) - t104, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * qJD(4) - t106, t135; 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t62 * qJD(2) + t101, 0, 0, 0, 0, 0, t52, -t49, 0, 0, 0, 0, 0, t81, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, -t100, 0, 0, 0, 0, 0, t149, -t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t111, 0, 0, 0, -t84 * t58, t84 * t56, t141, -t150, 0, -t109, 0, -t33 * qJD(3) + t117, t116 + t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * qJD(5), 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t100, 0, 0, 0, 0, 0, -t149, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76 * t121, -t77 * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, 0, -t76 * t85 + t82, -t77 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, -t150, 0, t32 * qJD(4), 0, -t15 * qJD(3) - t91, t92 + t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 * qJD(4), 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, 0, t76 * t120 - t82, t77 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
