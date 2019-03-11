% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:18
% EndTime: 2019-03-09 02:06:18
% DurationCPUTime: 0.10s
% Computational Cost: add. (227->41), mult. (400->84), div. (0->0), fcn. (191->6), ass. (0->32)
t144 = qJD(1) ^ 2;
t154 = t144 / 0.2e1;
t153 = sin(qJ(5));
t131 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t138 = sin(pkin(9));
t139 = cos(pkin(9));
t150 = qJ(2) * qJD(1);
t126 = t138 * t131 + t139 * t150;
t124 = -qJD(1) * pkin(7) + t126;
t141 = sin(qJ(4));
t143 = cos(qJ(4));
t151 = t141 * qJD(3) + t143 * t124;
t119 = qJD(4) * pkin(8) + t151;
t125 = t139 * t131 - t138 * t150;
t123 = qJD(1) * pkin(3) - t125;
t120 = (pkin(4) * t143 + pkin(8) * t141) * qJD(1) + t123;
t142 = cos(qJ(5));
t152 = t142 * t119 + t153 * t120;
t149 = qJD(1) * t141;
t148 = t143 * qJD(1);
t147 = qJD(1) * qJD(4);
t146 = t143 * qJD(3) - t141 * t124;
t118 = -qJD(4) * pkin(4) - t146;
t145 = -t119 * t153 + t142 * t120;
t135 = -qJD(1) * pkin(1) + qJD(2);
t132 = qJD(5) + t148;
t128 = -qJD(4) * t153 + t142 * t149;
t127 = t142 * qJD(4) + t149 * t153;
t115 = -t127 * pkin(5) + t128 * qJ(6) + t118;
t114 = t132 * qJ(6) + t152;
t113 = -t132 * pkin(5) + qJD(6) - t145;
t1 = [t154, 0, 0, -t135 * qJD(1), t144 * qJ(2), qJ(2) ^ 2 * t154 + t135 ^ 2 / 0.2e1, -t125 * qJD(1), t126 * qJD(1), t126 ^ 2 / 0.2e1 + t125 ^ 2 / 0.2e1 + qJD(3) ^ 2 / 0.2e1, t141 ^ 2 * t154, t141 * t144 * t143, -t141 * t147, -t143 * t147, qJD(4) ^ 2 / 0.2e1, qJD(4) * t146 + t123 * t148, -qJD(4) * t151 - t123 * t149, t128 ^ 2 / 0.2e1, -t128 * t127, -t128 * t132, t127 * t132, t132 ^ 2 / 0.2e1, -t118 * t127 + t132 * t145, -t118 * t128 - t132 * t152, -t113 * t132 - t115 * t127, -t113 * t128 + t114 * t127, t114 * t132 + t115 * t128, t114 ^ 2 / 0.2e1 + t115 ^ 2 / 0.2e1 + t113 ^ 2 / 0.2e1;];
T_reg  = t1;
