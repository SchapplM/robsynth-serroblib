% Calculate minimal parameter regressor of fixed base kinetic energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S4RRRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energykin_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_energykin_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:50
% EndTime: 2019-12-31 17:30:50
% DurationCPUTime: 0.07s
% Computational Cost: add. (145->31), mult. (385->73), div. (0->0), fcn. (283->8), ass. (0->33)
t135 = sin(pkin(4));
t143 = qJD(1) ^ 2;
t154 = t135 ^ 2 * t143;
t150 = cos(pkin(4)) * qJD(1);
t133 = qJD(2) + t150;
t139 = sin(qJ(2));
t142 = cos(qJ(2));
t151 = qJD(1) * t135;
t146 = t142 * t151;
t149 = pkin(1) * t150;
t152 = pkin(6) * t146 + t139 * t149;
t121 = t133 * pkin(7) + t152;
t122 = (-pkin(2) * t142 - pkin(7) * t139 - pkin(1)) * t151;
t138 = sin(qJ(3));
t141 = cos(qJ(3));
t153 = t141 * t121 + t138 * t122;
t148 = t142 * t154;
t147 = t139 * t151;
t145 = -t138 * t121 + t141 * t122;
t144 = -pkin(6) * t147 + t142 * t149;
t124 = -t141 * t133 + t138 * t147;
t120 = -t133 * pkin(2) - t144;
t140 = cos(qJ(4));
t137 = sin(qJ(4));
t128 = -qJD(3) + t146;
t125 = t138 * t133 + t141 * t147;
t123 = qJD(4) + t124;
t117 = t140 * t125 - t137 * t128;
t116 = t137 * t125 + t140 * t128;
t115 = -t128 * pkin(8) + t153;
t114 = t128 * pkin(3) - t145;
t113 = t124 * pkin(3) - t125 * pkin(8) + t120;
t1 = [t143 / 0.2e1, 0, 0, t139 ^ 2 * t154 / 0.2e1, t139 * t148, t133 * t147, t133 * t146, t133 ^ 2 / 0.2e1, pkin(1) * t148 + t144 * t133, -pkin(1) * t139 * t154 - t152 * t133, t125 ^ 2 / 0.2e1, -t125 * t124, -t125 * t128, t124 * t128, t128 ^ 2 / 0.2e1, t120 * t124 - t145 * t128, t120 * t125 + t153 * t128, t117 ^ 2 / 0.2e1, -t117 * t116, t117 * t123, -t116 * t123, t123 ^ 2 / 0.2e1, (t140 * t113 - t137 * t115) * t123 + t114 * t116, -(t137 * t113 + t140 * t115) * t123 + t114 * t117;];
T_reg = t1;
