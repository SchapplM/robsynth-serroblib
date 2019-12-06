% Calculate minimal parameter regressor of coriolis matrix for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5PRRRP3_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:15
% EndTime: 2019-12-05 16:44:18
% DurationCPUTime: 0.88s
% Computational Cost: add. (1436->100), mult. (2870->146), div. (0->0), fcn. (2978->4), ass. (0->90)
t155 = pkin(6) + pkin(7);
t93 = sin(qJ(3));
t85 = t155 * t93;
t92 = sin(qJ(4));
t146 = t92 * t85;
t149 = cos(qJ(4));
t94 = cos(qJ(3));
t86 = t155 * t94;
t84 = t149 * t86;
t160 = -t84 + t146;
t79 = -t149 * t94 + t92 * t93;
t40 = t79 * qJ(5) + t160;
t151 = t40 * pkin(4);
t113 = t149 * pkin(3);
t89 = t113 + pkin(4);
t99 = t113 / 0.2e1 - t89 / 0.2e1;
t119 = qJD(3) + qJD(4);
t46 = t149 * t85 + t92 * t86;
t81 = t149 * t93 + t92 * t94;
t162 = -t81 * qJ(5) - t46;
t161 = t119 * t79;
t159 = t119 * t46;
t78 = t79 ^ 2;
t158 = t81 ^ 2;
t102 = -t84 / 0.2e1;
t154 = pkin(3) * t92;
t153 = pkin(4) * t79;
t152 = pkin(4) * t81;
t150 = t93 * pkin(3);
t148 = t89 * t81;
t147 = t92 * t79;
t139 = qJD(2) * t81;
t90 = -t94 * pkin(3) - pkin(2);
t138 = qJD(2) * t90;
t137 = qJD(2) * t94;
t11 = -t162 * t81 + t79 * t40;
t136 = t11 * qJD(2);
t108 = -t147 / 0.2e1;
t109 = -t148 / 0.2e1;
t69 = -t152 / 0.2e1;
t18 = t109 + t69 + (t108 - t93 / 0.2e1) * pkin(3);
t133 = t18 * qJD(2);
t29 = (-pkin(4) / 0.2e1 - t99) * t79;
t132 = t29 * qJD(2);
t45 = t78 - t158;
t131 = t45 * qJD(2);
t48 = t79 * t150 + t90 * t81;
t128 = t48 * qJD(2);
t49 = t81 * t150 - t90 * t79;
t127 = t49 * qJD(2);
t53 = t78 + t158;
t126 = t53 * qJD(2);
t60 = t102 + t84 / 0.2e1;
t125 = t60 * qJD(2);
t124 = t79 * qJD(4);
t123 = t81 * qJD(4);
t87 = -t93 ^ 2 + t94 ^ 2;
t122 = t87 * qJD(2);
t121 = t93 * qJD(3);
t120 = t94 * qJD(3);
t118 = pkin(2) * t93 * qJD(2);
t117 = pkin(2) * t137;
t116 = pkin(4) * t139;
t115 = pkin(4) * t123;
t114 = qJD(4) * t154;
t112 = t79 * t138;
t111 = t81 * t138;
t110 = t93 * t137;
t105 = t149 * qJD(3);
t104 = t149 * qJD(4);
t52 = t119 * t81;
t63 = t90 + t153;
t9 = t63 * t152;
t101 = t9 * qJD(2);
t8 = t63 * (t150 + t152);
t100 = t8 * qJD(2);
t97 = t99 * t81;
t95 = t99 * t40;
t3 = -t151 / 0.2e1 - t95;
t70 = t152 / 0.2e1;
t31 = t70 + t97;
t67 = (t113 - t89) * t154;
t96 = -t31 * qJD(1) - t3 * qJD(2) - t67 * qJD(3);
t58 = t79 * t139;
t47 = 0.2e1 * t102 + t146;
t30 = t69 + t97;
t28 = t153 / 0.2e1 - t99 * t79;
t17 = pkin(3) * t108 + t109 + t150 / 0.2e1 + t70;
t2 = t151 / 0.2e1 - t95;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t120, 0, 0, 0, 0, 0, -t52, t161, 0, (-pkin(3) * t147 - t148) * qJD(3) + t30 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t161, 0, t30 * qJD(3) - t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, t93 * t120, t87 * qJD(3), 0, 0, 0, -pkin(2) * t121, -pkin(2) * t120, -t79 * t52, t119 * t45, 0, 0, 0, t48 * qJD(3) + t90 * t123, t49 * qJD(3) - t90 * t124, t53 * qJD(5), t8 * qJD(3) + t9 * qJD(4) + t11 * qJD(5); 0, 0, 0, 0, t110, t122, t120, -t121, 0, -pkin(6) * t120 - t118, pkin(6) * t121 - t117, -t58, t131, -t161, -t52, 0, qJD(3) * t160 + t47 * qJD(4) + t128, t127 + t159, (-t81 * t154 + t89 * t79) * qJD(3) + t28 * qJD(4), (t154 * t162 + t40 * t89) * qJD(3) + t2 * qJD(4) + t17 * qJD(5) + t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t131, -t161, -t52, 0, t47 * qJD(3) + qJD(4) * t160 + t111, -t112 + t159, pkin(4) * t124 + t28 * qJD(3), t2 * qJD(3) + qJD(4) * t151 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, t17 * qJD(3) + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * qJD(4); 0, 0, 0, 0, -t110, -t122, 0, 0, 0, t118, t117, t58, -t131, 0, 0, 0, t60 * qJD(4) - t128, -t127, t29 * qJD(4), t3 * qJD(4) + t18 * qJD(5) - t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, -pkin(3) * t104, 0, t67 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119 * t154 + t125, (-t105 - t104) * pkin(3), t132, -pkin(4) * t114 - t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t131, 0, 0, 0, -t60 * qJD(3) - t111, t112, -t29 * qJD(3), -t3 * qJD(3) - qJD(5) * t152 - t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t154 - t125, pkin(3) * t105, -t132, t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, -t18 * qJD(3) + t115 - t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t1;
