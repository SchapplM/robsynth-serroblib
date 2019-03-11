% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:32:01
% EndTime: 2019-03-09 01:32:01
% DurationCPUTime: 0.08s
% Computational Cost: add. (163->40), mult. (351->85), div. (0->0), fcn. (203->8), ass. (0->33)
t139 = cos(pkin(9));
t147 = -pkin(1) * t139 - pkin(2);
t126 = qJD(3) + (-qJ(4) + t147) * qJD(1);
t136 = sin(pkin(10));
t138 = cos(pkin(10));
t118 = -t136 * qJD(2) + t138 * t126;
t148 = qJD(1) * t138;
t114 = -pkin(7) * t148 + t118;
t119 = t138 * qJD(2) + t136 * t126;
t149 = qJD(1) * t136;
t115 = -pkin(7) * t149 + t119;
t141 = sin(qJ(5));
t143 = cos(qJ(5));
t150 = t141 * t114 + t143 * t115;
t137 = sin(pkin(9));
t129 = (-pkin(1) * t137 - qJ(3)) * qJD(1);
t127 = qJD(4) - t129;
t123 = pkin(4) * t149 + t127;
t146 = t143 * t114 - t141 * t115;
t124 = (t136 * t143 + t138 * t141) * qJD(1);
t144 = qJD(1) ^ 2;
t142 = cos(qJ(6));
t140 = sin(qJ(6));
t134 = qJD(2) ^ 2 / 0.2e1;
t128 = t147 * qJD(1) + qJD(3);
t125 = (-t136 * t141 + t138 * t143) * qJD(1);
t122 = qJD(6) + t124;
t117 = t140 * qJD(5) + t142 * t125;
t116 = -t142 * qJD(5) + t140 * t125;
t111 = t124 * pkin(5) - t125 * pkin(8) + t123;
t110 = qJD(5) * pkin(8) + t150;
t109 = -qJD(5) * pkin(5) - t146;
t1 = [t144 / 0.2e1, 0, 0, t134 + (t137 ^ 2 / 0.2e1 + t139 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t144, t128 * qJD(1), -t129 * qJD(1), t134 + t129 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1, t127 * t149, t127 * t148 (-t118 * t138 - t119 * t136) * qJD(1), t119 ^ 2 / 0.2e1 + t118 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1, t125 ^ 2 / 0.2e1, -t125 * t124, t125 * qJD(5), -t124 * qJD(5), qJD(5) ^ 2 / 0.2e1, t146 * qJD(5) + t123 * t124, -t150 * qJD(5) + t123 * t125, t117 ^ 2 / 0.2e1, -t117 * t116, t117 * t122, -t116 * t122, t122 ^ 2 / 0.2e1 (-t140 * t110 + t142 * t111) * t122 + t109 * t116 -(t142 * t110 + t140 * t111) * t122 + t109 * t117;];
T_reg  = t1;
