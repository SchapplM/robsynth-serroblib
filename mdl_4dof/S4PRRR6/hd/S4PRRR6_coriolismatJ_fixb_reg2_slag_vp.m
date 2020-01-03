% Calculate inertial parameters regressor of coriolis matrix for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4PRRR6_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:09
% EndTime: 2019-12-31 16:35:10
% DurationCPUTime: 0.80s
% Computational Cost: add. (1096->85), mult. (2806->148), div. (0->0), fcn. (2883->6), ass. (0->81)
t126 = qJD(3) + qJD(4);
t104 = sin(qJ(4));
t105 = sin(qJ(3));
t155 = cos(qJ(4));
t120 = t155 * t105;
t157 = -pkin(6) - pkin(5);
t107 = cos(qJ(3));
t92 = t157 * t107;
t115 = -t104 * t92 - t157 * t120;
t108 = cos(qJ(2));
t119 = t155 * t107;
t144 = t104 * t105;
t109 = -t119 / 0.2e1 + t144 / 0.2e1;
t112 = t119 - t144;
t47 = (t112 / 0.2e1 + t109) * t108;
t136 = t47 * qJD(1);
t162 = t126 * t115 - t136;
t66 = t157 * t144 - t155 * t92;
t161 = t126 * t66;
t143 = t104 * t107;
t86 = t120 + t143;
t101 = -t107 * pkin(3) - pkin(2);
t133 = qJD(2) * t101;
t110 = -t143 / 0.2e1 - t120 / 0.2e1;
t46 = (t86 / 0.2e1 + t110) * t108;
t137 = t46 * qJD(1);
t160 = -t86 * t133 + t137;
t153 = pkin(3) * t105;
t44 = t101 * t86 - t112 * t153;
t159 = -t44 * qJD(2) + t137;
t154 = pkin(3) * t104;
t150 = t112 * t86;
t148 = qJD(3) * pkin(3);
t142 = t108 * t105;
t106 = sin(qJ(2));
t74 = t86 * t106;
t75 = t108 * t86;
t76 = t112 * t106;
t77 = t112 * t108;
t99 = t106 * t108;
t16 = t74 * t75 + t76 * t77 - t99;
t141 = t16 * qJD(1);
t34 = t112 ^ 2 - t86 ^ 2;
t140 = t34 * qJD(2);
t45 = t101 * t112 + t86 * t153;
t138 = t45 * qJD(2);
t102 = t105 ^ 2;
t103 = t107 ^ 2;
t114 = (t102 + t103) * t108;
t67 = t106 * t114 - t99;
t135 = t67 * qJD(1);
t97 = t103 - t102;
t134 = t97 * qJD(2);
t132 = qJD(2) * t107;
t131 = qJD(4) * t101;
t130 = t105 * qJD(3);
t129 = t106 * qJD(2);
t128 = t107 * qJD(3);
t127 = t108 * qJD(2);
t125 = pkin(2) * t105 * qJD(2);
t124 = pkin(2) * t132;
t123 = qJD(2) * t150;
t122 = t112 * t133;
t118 = t105 * t128;
t117 = t155 * qJD(3);
t116 = t155 * qJD(4);
t58 = t126 * t86;
t111 = t77 * t104 / 0.2e1 - t75 * t155 / 0.2e1;
t1 = (t142 / 0.2e1 + t111) * pkin(3);
t14 = t101 * t153;
t113 = -t1 * qJD(1) + t14 * qJD(2);
t98 = t105 * t132;
t57 = t126 * t112;
t49 = (-t86 / 0.2e1 + t110) * t108;
t48 = (-t112 / 0.2e1 + t109) * t108;
t39 = t47 * qJD(2);
t37 = t46 * qJD(2);
t18 = t49 * qJD(2) - t126 * t76;
t17 = t48 * qJD(2) + t126 * t74;
t2 = (-t142 / 0.2e1 + t111) * pkin(3);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, -t127, 0, 0, 0, 0, 0, 0, 0, 0, -t107 * t129 - t108 * t130, t105 * t129 - t108 * t128, qJD(2) * t114, t135 + (-t106 * pkin(2) + pkin(5) * t114) * qJD(2), 0, 0, 0, 0, 0, 0, -t112 * t129 + t126 * t49, t126 * t48 + t129 * t86, (t112 * t77 + t75 * t86) * qJD(2), t141 + (t106 * t101 + t115 * t75 + t77 * t66) * qJD(2) + t2 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105 * t127 - t106 * t128, t106 * t130 - t107 * t127, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17, 0, t2 * qJD(2) + (-t104 * t74 - t155 * t76) * t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, 0, 0, 0, 0, 0, 0, -t126 * t46, -t126 * t47, 0, -t1 * qJD(3) - t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t97 * qJD(3), 0, -t118, 0, 0, -pkin(2) * t130, -pkin(2) * t128, 0, 0, t112 * t58, t126 * t34, 0, -t126 * t150, 0, 0, t44 * qJD(3) + t131 * t86, t45 * qJD(3) + t112 * t131, 0, t14 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t134, t128, -t98, -t130, 0, -pkin(5) * t128 - t125, pkin(5) * t130 - t124, 0, 0, t123, t140, t57, -t123, -t58, 0, -t159 - t161, t138 + t162, (-t104 * t86 - t112 * t155) * t148, (-t104 * t115 - t155 * t66) * t148 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, t140, t57, -t123, -t58, 0, -t160 - t161, t122 + t162, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t39, 0, t1 * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t134, 0, t98, 0, 0, t125, t124, 0, 0, -t123, -t140, 0, t123, 0, 0, t159, t136 - t138, 0, -t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t154, -pkin(3) * t116, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 * t154, (-t117 - t116) * pkin(3), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, -t140, 0, t123, 0, 0, t160, t136 - t122, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104 * t148, pkin(3) * t117, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
