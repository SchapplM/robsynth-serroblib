% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:47
% EndTime: 2019-03-09 02:03:47
% DurationCPUTime: 0.07s
% Computational Cost: add. (174->39), mult. (335->80), div. (0->0), fcn. (164->6), ass. (0->32)
t133 = qJD(1) ^ 2;
t145 = t133 / 0.2e1;
t128 = cos(pkin(9));
t138 = -pkin(1) * t128 - pkin(2);
t116 = qJD(3) + (-pkin(7) + t138) * qJD(1);
t130 = sin(qJ(4));
t132 = cos(qJ(4));
t143 = t132 * qJD(2) + t130 * t116;
t112 = qJD(4) * pkin(8) + t143;
t127 = sin(pkin(9));
t137 = -pkin(1) * t127 - qJ(3);
t114 = (pkin(4) * t130 - pkin(8) * t132 - t137) * qJD(1);
t129 = sin(qJ(5));
t131 = cos(qJ(5));
t144 = t131 * t112 + t129 * t114;
t142 = qJD(1) * t132;
t120 = t137 * qJD(1);
t141 = t120 * qJD(1);
t140 = t130 * qJD(1);
t139 = qJD(1) * qJD(4);
t136 = -t130 * qJD(2) + t132 * t116;
t135 = -t129 * t112 + t131 * t114;
t111 = -qJD(4) * pkin(4) - t136;
t126 = qJD(2) ^ 2 / 0.2e1;
t122 = qJD(5) + t140;
t119 = t138 * qJD(1) + qJD(3);
t118 = t129 * qJD(4) + t131 * t142;
t117 = -t131 * qJD(4) + t129 * t142;
t109 = t117 * pkin(5) - t118 * qJ(6) + t111;
t108 = t122 * qJ(6) + t144;
t107 = -t122 * pkin(5) + qJD(6) - t135;
t1 = [t145, 0, 0, t126 + (t127 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t133, t119 * qJD(1), -t141, t126 + t120 ^ 2 / 0.2e1 + t119 ^ 2 / 0.2e1, t132 ^ 2 * t145, -t132 * t133 * t130, t132 * t139, -t130 * t139, qJD(4) ^ 2 / 0.2e1, t136 * qJD(4) - t120 * t140, -t143 * qJD(4) - t132 * t141, t118 ^ 2 / 0.2e1, -t118 * t117, t118 * t122, -t117 * t122, t122 ^ 2 / 0.2e1, t111 * t117 + t135 * t122, t111 * t118 - t144 * t122, -t107 * t122 + t109 * t117, t107 * t118 - t108 * t117, t108 * t122 - t109 * t118, t108 ^ 2 / 0.2e1 + t109 ^ 2 / 0.2e1 + t107 ^ 2 / 0.2e1;];
T_reg  = t1;
