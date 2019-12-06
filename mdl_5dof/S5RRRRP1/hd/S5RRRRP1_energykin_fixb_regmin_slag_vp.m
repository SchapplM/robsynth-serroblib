% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:45:59
% EndTime: 2019-12-05 18:45:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (194->35), mult. (494->78), div. (0->0), fcn. (338->6), ass. (0->35)
t138 = -pkin(7) - pkin(6);
t125 = qJD(1) ^ 2;
t137 = t125 / 0.2e1;
t136 = cos(qJ(4));
t124 = cos(qJ(2));
t135 = t124 * t125;
t121 = sin(qJ(3));
t122 = sin(qJ(2));
t123 = cos(qJ(3));
t112 = (t121 * t124 + t122 * t123) * qJD(1);
t119 = qJD(2) + qJD(3);
t132 = qJD(1) * t122;
t114 = qJD(2) * pkin(2) + t138 * t132;
t131 = qJD(1) * t124;
t115 = t138 * t131;
t126 = t123 * t114 + t115 * t121;
t102 = pkin(3) * t119 - pkin(8) * t112 + t126;
t111 = t121 * t132 - t123 * t131;
t133 = t121 * t114 - t123 * t115;
t104 = -pkin(8) * t111 + t133;
t120 = sin(qJ(4));
t134 = t120 * t102 + t136 * t104;
t130 = qJD(1) * qJD(2);
t129 = t122 * t130;
t128 = t124 * t130;
t127 = t136 * t102 - t104 * t120;
t116 = (-pkin(2) * t124 - pkin(1)) * qJD(1);
t107 = pkin(3) * t111 + t116;
t118 = qJD(4) + t119;
t106 = -t120 * t111 + t136 * t112;
t105 = t136 * t111 + t112 * t120;
t99 = pkin(4) * t105 + qJD(5) + t107;
t98 = -qJ(5) * t105 + t134;
t97 = pkin(4) * t118 - qJ(5) * t106 + t127;
t1 = [t137, 0, 0, t122 ^ 2 * t137, t122 * t135, t129, t128, qJD(2) ^ 2 / 0.2e1, pkin(1) * t135 - pkin(6) * t129, -pkin(1) * t122 * t125 - pkin(6) * t128, t112 ^ 2 / 0.2e1, -t112 * t111, t112 * t119, -t111 * t119, t119 ^ 2 / 0.2e1, t116 * t111 + t126 * t119, t116 * t112 - t133 * t119, t106 ^ 2 / 0.2e1, -t106 * t105, t106 * t118, -t105 * t118, t118 ^ 2 / 0.2e1, t107 * t105 + t127 * t118, t107 * t106 - t134 * t118, -t105 * t98 - t106 * t97, t98 ^ 2 / 0.2e1 + t97 ^ 2 / 0.2e1 + t99 ^ 2 / 0.2e1;];
T_reg = t1;
