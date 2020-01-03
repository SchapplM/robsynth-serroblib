% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRP8
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
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:26
% EndTime: 2019-12-31 20:04:26
% DurationCPUTime: 0.07s
% Computational Cost: add. (107->34), mult. (257->75), div. (0->0), fcn. (126->4), ass. (0->28)
t127 = qJD(1) ^ 2;
t137 = t127 / 0.2e1;
t126 = cos(qJ(2));
t136 = t126 * t127;
t124 = sin(qJ(2));
t134 = qJD(1) * t124;
t132 = pkin(6) * t134 + qJD(3);
t106 = -pkin(7) * t134 + (-pkin(2) - pkin(3)) * qJD(2) + t132;
t133 = qJD(1) * t126;
t113 = pkin(6) * t133 + qJD(2) * qJ(3);
t110 = -pkin(7) * t133 + t113;
t123 = sin(qJ(4));
t125 = cos(qJ(4));
t135 = t123 * t106 + t125 * t110;
t131 = qJD(1) * qJD(2);
t111 = -qJD(1) * pkin(1) - pkin(2) * t133 - qJ(3) * t134;
t130 = t124 * t131;
t129 = t126 * t131;
t128 = t125 * t106 - t123 * t110;
t105 = pkin(3) * t133 - t111;
t119 = qJD(2) - qJD(4);
t112 = -qJD(2) * pkin(2) + t132;
t109 = (-t123 * t126 + t124 * t125) * qJD(1);
t108 = (t123 * t124 + t125 * t126) * qJD(1);
t102 = t108 * pkin(4) + qJD(5) + t105;
t101 = -t108 * qJ(5) + t135;
t100 = -t119 * pkin(4) - t109 * qJ(5) + t128;
t1 = [t137, 0, 0, t124 ^ 2 * t137, t124 * t136, t130, t129, qJD(2) ^ 2 / 0.2e1, pkin(1) * t136 - pkin(6) * t130, -t127 * pkin(1) * t124 - pkin(6) * t129, -t112 * qJD(2) - t111 * t133, (t112 * t124 + t113 * t126) * qJD(1), t113 * qJD(2) - t111 * t134, t113 ^ 2 / 0.2e1 + t111 ^ 2 / 0.2e1 + t112 ^ 2 / 0.2e1, t109 ^ 2 / 0.2e1, -t109 * t108, -t109 * t119, t108 * t119, t119 ^ 2 / 0.2e1, t105 * t108 - t128 * t119, t105 * t109 + t135 * t119, -t100 * t109 - t101 * t108, t101 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1 + t102 ^ 2 / 0.2e1;];
T_reg = t1;
