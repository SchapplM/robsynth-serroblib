% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x25]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR10_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:20
% EndTime: 2019-12-31 18:04:23
% DurationCPUTime: 0.94s
% Computational Cost: add. (919->101), mult. (1885->149), div. (0->0), fcn. (2152->6), ass. (0->93)
t104 = qJD(4) + qJD(5);
t80 = sin(pkin(8));
t81 = cos(pkin(8));
t83 = sin(qJ(4));
t85 = cos(qJ(4));
t58 = t80 * t83 + t81 * t85;
t133 = -pkin(6) + qJ(2);
t68 = t133 * t80;
t69 = t133 * t81;
t87 = -t83 * t68 - t69 * t85;
t29 = -t58 * pkin(7) - t87;
t82 = sin(qJ(5));
t84 = cos(qJ(5));
t60 = t80 * t85 - t81 * t83;
t92 = -t85 * t68 + t69 * t83;
t86 = -pkin(7) * t60 - t92;
t159 = t104 * (-t84 * t29 - t82 * t86);
t158 = t104 * (t82 * t29 - t84 * t86);
t138 = t82 * t58;
t54 = t84 * t60;
t146 = t54 - t138;
t137 = t82 * t60;
t53 = t84 * t58;
t36 = t53 + t137;
t150 = -t146 ^ 2 + t36 ^ 2;
t157 = t150 * qJD(1);
t144 = t146 * qJD(1);
t154 = t36 * t144;
t143 = -t54 / 0.2e1;
t50 = t54 / 0.2e1;
t147 = t143 + t50;
t151 = qJD(4) * t147;
t145 = qJD(5) * t146;
t48 = -t53 / 0.2e1;
t49 = t53 / 0.2e1;
t142 = pkin(4) * t60;
t33 = t49 + t48;
t78 = t80 ^ 2;
t73 = t81 ^ 2 + t78;
t132 = pkin(4) * qJD(5);
t129 = qJD(4) * pkin(4);
t66 = -t81 * pkin(2) - t80 * qJ(3) - pkin(1);
t55 = t81 * pkin(3) - t66;
t40 = pkin(4) * t58 + t55;
t6 = -t36 * t142 - t146 * t40;
t127 = t6 * qJD(1);
t7 = -t142 * t146 + t36 * t40;
t126 = t7 * qJD(1);
t9 = -t138 + 0.2e1 * t50;
t124 = t9 * qJD(1);
t122 = qJD(1) * t40;
t121 = qJD(3) * t80;
t120 = qJD(4) * t36;
t119 = qJD(5) * t36;
t118 = qJD(5) * t40;
t46 = -t137 / 0.2e1;
t90 = 0.2e1 * t46;
t11 = 0.2e1 * t48 + t90;
t117 = t11 * qJD(1);
t17 = t137 + 0.2e1 * t49;
t115 = t17 * qJD(1);
t30 = t58 ^ 2 - t60 ^ 2;
t114 = t30 * qJD(1);
t32 = t46 + t137 / 0.2e1;
t113 = t32 * qJD(1);
t112 = t147 * qJD(1);
t111 = t58 * qJD(1);
t110 = t58 * qJD(4);
t109 = t60 * qJD(1);
t56 = t60 * qJD(4);
t67 = t73 * qJ(2);
t108 = t67 * qJD(1);
t107 = t73 * qJD(1);
t106 = t78 * qJD(1);
t105 = t80 * qJD(1);
t101 = t36 * t122;
t100 = t146 * t122;
t99 = t58 * t109;
t98 = t36 * t105;
t97 = t146 * t105;
t96 = t58 * t105;
t95 = t60 * t105;
t94 = t81 * t105;
t91 = pkin(4) * t104;
t89 = qJD(1) * t55 - qJD(2);
t70 = t73 * qJD(2);
t63 = t67 * qJD(2);
t42 = t104 * (t82 * t83 - t84 * t85);
t41 = t104 * (-t82 * t85 - t84 * t83);
t18 = 0.2e1 * t143 + t138;
t16 = -t53 + t90;
t12 = t32 + t33;
t1 = [0, 0, 0, 0, 0, t70, t63, t81 * t121, t70, t78 * qJD(3), -t66 * t121 + t63, -t58 * t56, t30 * qJD(4), 0, 0, 0, t58 * t121 + t55 * t56, -t55 * t110 + t60 * t121, (-t119 - t120) * t146, t104 * t150, 0, 0, 0, -qJD(4) * t6 + t118 * t146 + t36 * t121, -qJD(4) * t7 - t36 * t118 + t121 * t146; 0, 0, 0, 0, 0, t107, t108, 0, t107, 0, t108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, qJD(4) * t12 + qJD(5) * t33; 0, 0, 0, 0, 0, 0, 0, t94, 0, t106, -t66 * t105, 0, 0, 0, 0, 0, t96, t95, 0, 0, 0, 0, 0, t98, t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t114, -t110, -t56, 0, t87 * qJD(4) + t55 * t109, t92 * qJD(4) - t55 * t111, -t154, t157, qJD(5) * t16 - t120, -qJD(4) * t146 + qJD(5) * t18, 0, t147 * qJD(2) - t127 + t159, t12 * qJD(2) - t126 + t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t157, qJD(4) * t16 - t119, qJD(4) * t18 - t145, 0, t100 + t159, t33 * qJD(2) - t101 + t158; 0, 0, 0, 0, 0, -t107, -t108, 0, -t107, 0, -t108 - t121, 0, 0, 0, 0, 0, -t56, t110, 0, 0, 0, 0, 0, -qJD(4) * t9 - t145, -qJD(4) * t11 + qJD(5) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, t111, 0, 0, 0, 0, 0, -t124, -t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t115; 0, 0, 0, 0, 0, 0, 0, -t94, 0, -t106, (qJD(1) * t66 + qJD(2)) * t80, 0, 0, 0, 0, 0, -t96, -t95, 0, 0, 0, 0, 0, -t98, -t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83 * qJD(4), -t85 * qJD(4), 0, 0, 0, 0, 0, t41, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, -t114, 0, 0, 0, -t89 * t60, t89 * t58, t154, -t157, t32 * qJD(5), t147 * qJD(5), 0, qJD(2) * t9 + t127, qJD(2) * t11 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, -t111, 0, 0, 0, 0, 0, t124, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82 * t132, -t84 * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t112, 0, -t82 * t91, -t84 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, -t157, -t32 * qJD(4), -t151, 0, qJD(2) * t146 - t100, -qJD(2) * t17 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, -t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, -t112, 0, t82 * t129, t84 * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
