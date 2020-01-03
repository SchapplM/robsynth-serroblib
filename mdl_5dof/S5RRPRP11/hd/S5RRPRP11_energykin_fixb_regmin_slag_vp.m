% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:03
% EndTime: 2019-12-31 20:14:03
% DurationCPUTime: 0.12s
% Computational Cost: add. (148->36), mult. (320->78), div. (0->0), fcn. (155->4), ass. (0->30)
t140 = -pkin(2) - pkin(7);
t128 = qJD(1) ^ 2;
t139 = t128 / 0.2e1;
t127 = cos(qJ(2));
t138 = t127 * t128;
t125 = sin(qJ(2));
t130 = -qJ(3) * t125 - pkin(1);
t111 = (t140 * t127 + t130) * qJD(1);
t135 = t125 * qJD(1);
t134 = pkin(6) * t135 + qJD(3);
t112 = pkin(3) * t135 + t140 * qJD(2) + t134;
t124 = sin(qJ(4));
t126 = cos(qJ(4));
t137 = t126 * t111 + t124 * t112;
t136 = qJD(1) * t127;
t118 = -pkin(6) * t136 - qJD(2) * qJ(3);
t133 = qJD(1) * qJD(2);
t113 = pkin(3) * t136 - t118;
t132 = t125 * t133;
t131 = t127 * t133;
t129 = -t124 * t111 + t126 * t112;
t119 = qJD(4) + t135;
t117 = -qJD(2) * pkin(2) + t134;
t116 = t126 * qJD(2) - t124 * t136;
t115 = t124 * qJD(2) + t126 * t136;
t114 = (-pkin(2) * t127 + t130) * qJD(1);
t108 = t115 * pkin(4) - t116 * qJ(5) + t113;
t107 = t119 * qJ(5) + t137;
t106 = -t119 * pkin(4) + qJD(5) - t129;
t1 = [t139, 0, 0, t125 ^ 2 * t139, t125 * t138, t132, t131, qJD(2) ^ 2 / 0.2e1, pkin(1) * t138 - pkin(6) * t132, -t128 * pkin(1) * t125 - pkin(6) * t131, (t117 * t125 - t118 * t127) * qJD(1), t117 * qJD(2) + t114 * t136, -t118 * qJD(2) - t114 * t135, t114 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1 + t117 ^ 2 / 0.2e1, t116 ^ 2 / 0.2e1, -t116 * t115, t116 * t119, -t115 * t119, t119 ^ 2 / 0.2e1, t113 * t115 + t129 * t119, t113 * t116 - t137 * t119, -t106 * t119 + t108 * t115, t106 * t116 - t107 * t115, t107 * t119 - t108 * t116, t107 ^ 2 / 0.2e1 + t108 ^ 2 / 0.2e1 + t106 ^ 2 / 0.2e1;];
T_reg = t1;
