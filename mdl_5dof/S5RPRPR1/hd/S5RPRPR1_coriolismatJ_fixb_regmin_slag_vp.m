% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:34
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:33:57
% EndTime: 2021-01-15 11:34:02
% DurationCPUTime: 1.36s
% Computational Cost: add. (1532->121), mult. (2693->162), div. (0->0), fcn. (3035->6), ass. (0->104)
t203 = qJD(3) + qJD(5);
t105 = sin(qJ(5));
t107 = cos(qJ(5));
t103 = sin(pkin(8));
t104 = cos(pkin(8));
t106 = sin(qJ(3));
t109 = -pkin(1) - pkin(6);
t181 = -qJ(4) + t109;
t95 = t181 * t106;
t108 = cos(qJ(3));
t96 = t181 * t108;
t118 = -t103 * t95 + t104 * t96;
t91 = -t103 * t106 + t104 * t108;
t112 = -pkin(7) * t91 + t118;
t55 = -t103 * t96 - t104 * t95;
t92 = t103 * t108 + t104 * t106;
t38 = -pkin(7) * t92 - t55;
t202 = t203 * (-t105 * t112 - t107 * t38);
t201 = t203 * (t105 * t38 - t107 * t112);
t170 = t105 * t91;
t80 = t107 * t92;
t49 = t80 + t170;
t169 = t105 * t92;
t81 = t107 * t91;
t51 = t81 - t169;
t191 = t49 ^ 2 - t51 ^ 2;
t196 = t191 * qJD(1);
t193 = t51 * qJD(1);
t192 = t51 * qJD(3);
t142 = t49 * qJD(3);
t147 = t49 * qJD(5);
t15 = -t142 - t147;
t188 = qJD(4) * t49;
t187 = t51 * qJD(5);
t186 = t49 * qJD(1);
t185 = -t187 - t192;
t176 = qJD(3) * pkin(3);
t182 = (t103 * t91 - t104 * t92) * t176;
t180 = t91 ^ 2;
t179 = t92 ^ 2;
t123 = -t80 / 0.2e1;
t124 = t81 / 0.2e1;
t178 = pkin(3) * t103;
t177 = t108 * pkin(3);
t99 = pkin(3) * t106 + qJ(2);
t65 = pkin(4) * t92 + t99;
t66 = pkin(4) * t91 + t177;
t12 = t49 * t66 + t51 * t65;
t161 = t12 * qJD(1);
t13 = -t49 * t65 + t51 * t66;
t160 = t13 * qJD(1);
t16 = -t118 * t91 + t55 * t92;
t158 = t16 * qJD(1);
t117 = -t180 / 0.2e1 - t179 / 0.2e1;
t18 = -0.1e1 / 0.2e1 + t117;
t157 = t18 * qJD(1);
t22 = 0.2e1 * t124 - t169;
t154 = t22 * qJD(1);
t111 = -t103 * t92 / 0.2e1 - t104 * t91 / 0.2e1;
t41 = (-t108 / 0.2e1 + t111) * pkin(3);
t153 = t41 * qJD(1);
t42 = t177 * t91 - t92 * t99;
t152 = t42 * qJD(1);
t43 = t177 * t92 + t91 * t99;
t151 = t43 * qJD(1);
t47 = t124 - t81 / 0.2e1;
t150 = t47 * qJD(1);
t149 = t47 * qJD(5);
t58 = t179 + t180;
t141 = t58 * qJD(1);
t140 = t91 * qJD(1);
t139 = t91 * qJD(3);
t138 = t92 * qJD(1);
t97 = t106 ^ 2 - t108 ^ 2;
t137 = t97 * qJD(1);
t136 = t99 * qJD(1);
t135 = qJD(1) * qJ(2);
t134 = t106 * qJD(1);
t133 = t106 * qJD(3);
t132 = t108 * qJD(1);
t131 = t108 * qJD(3);
t130 = t49 * t193;
t129 = t51 * t186;
t128 = t65 * t186;
t127 = t65 * t193;
t121 = qJ(2) * t134;
t120 = qJ(2) * t132;
t119 = t106 * t132;
t9 = t99 * t177;
t116 = t9 * qJD(1);
t100 = pkin(3) * t104 + pkin(4);
t78 = -t100 * t107 + t105 * t178;
t114 = t78 * qJD(3);
t113 = qJD(5) * t22 + t192;
t48 = t123 + t80 / 0.2e1;
t79 = t100 * t105 + t107 * t178;
t110 = -qJD(2) * t48 + qJD(3) * t79;
t84 = t92 * qJD(3);
t68 = t79 * qJD(5);
t67 = t78 * qJD(5);
t40 = t177 / 0.2e1 + t111 * pkin(3);
t23 = 0.2e1 * t123 - t170;
t17 = 0.1e1 / 0.2e1 + t117;
t1 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t106 * t131, t97 * qJD(3), 0, 0, 0, qJ(2) * t131 + qJD(2) * t106, -qJ(2) * t133 + qJD(2) * t108, qJD(2) * t92 + qJD(3) * t43, qJD(2) * t91 + qJD(3) * t42, t58 * qJD(4), qJD(2) * t99 + qJD(3) * t9 + qJD(4) * t16, t15 * t51, t203 * t191, 0, 0, 0, qJD(2) * t49 + qJD(3) * t12 + t187 * t65, qJD(2) * t51 + qJD(3) * t13 - t147 * t65; 0, 0, 0, 0, qJD(1), t135, 0, 0, 0, 0, 0, t134, t132, t138, t140, 0, qJD(4) * t17 + t136, 0, 0, 0, 0, 0, t186, t193; 0, 0, 0, 0, 0, 0, -t119, t137, -t133, -t131, 0, -t109 * t133 + t120, -t109 * t131 - t121, qJD(3) * t55 + t151, -qJD(3) * t118 + t152, -t182, t40 * qJD(4) + (t103 * t118 + t104 * t55) * t176 + t116, -t129, t196, t15, -t113, 0, t161 + t202, t160 + t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, qJD(2) * t17 + qJD(3) * t40 + t158, 0, 0, 0, 0, 0, t149, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t196, t15, -qJD(3) * t22 - t187, 0, t47 * qJD(4) + t127 + t202, -t128 + t201; 0, 0, 0, 0, -qJD(1), -t135, 0, 0, 0, 0, 0, -t134, -t132, -t138, -t140, 0, qJD(4) * t18 - t136, 0, 0, 0, 0, 0, -t186, -t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t131, -t84, -t139, 0, t182, 0, 0, 0, 0, 0, qJD(5) * t23 - t142, t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t23 - t147, t185; 0, 0, 0, 0, 0, 0, t119, -t137, 0, 0, 0, -t120, t121, -qJD(4) * t91 - t151, qJD(4) * t92 - t152, 0, qJD(4) * t41 - t116, t129, -t196, 0, -t149, 0, -qJD(4) * t51 - t161, -t160 + t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * qJD(5), 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, t138, 0, t153, 0, 0, 0, 0, 0, -t193, t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, 0, -t110 - t68, t114 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, -t84, -t141, -qJD(2) * t18 - qJD(3) * t41 - t158, 0, 0, 0, 0, 0, t113, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t157, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t138, 0, -t153, 0, 0, 0, 0, 0, t193, -t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, -t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t196, 0, t47 * qJD(3), 0, -qJD(4) * t22 - t127, t128 + t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48 * qJD(3), 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, 0, t110, -t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
