% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:43
% EndTime: 2019-03-09 01:58:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (224->45), mult. (533->92), div. (0->0), fcn. (347->8), ass. (0->36)
t164 = cos(qJ(5));
t149 = sin(pkin(9));
t143 = (pkin(1) * t149 + qJ(3)) * qJD(1);
t150 = cos(pkin(10));
t146 = t150 * qJD(2);
t148 = sin(pkin(10));
t131 = t146 + (-pkin(7) * qJD(1) - t143) * t148;
t136 = t148 * qJD(2) + t150 * t143;
t160 = qJD(1) * t150;
t132 = pkin(7) * t160 + t136;
t153 = sin(qJ(4));
t154 = cos(qJ(4));
t162 = t153 * t131 + t154 * t132;
t124 = qJD(4) * pkin(8) + t162;
t151 = cos(pkin(9));
t159 = -pkin(1) * t151 - pkin(2);
t138 = qJD(3) + (-pkin(3) * t150 + t159) * qJD(1);
t161 = qJD(1) * t148;
t139 = t153 * t161 - t154 * t160;
t140 = (t148 * t154 + t150 * t153) * qJD(1);
t127 = t139 * pkin(4) - t140 * pkin(8) + t138;
t152 = sin(qJ(5));
t163 = t164 * t124 + t152 * t127;
t158 = -t152 * t124 + t164 * t127;
t157 = t154 * t131 - t153 * t132;
t123 = -qJD(4) * pkin(4) - t157;
t155 = qJD(1) ^ 2;
t142 = t159 * qJD(1) + qJD(3);
t137 = qJD(5) + t139;
t135 = -t148 * t143 + t146;
t134 = t152 * qJD(4) + t164 * t140;
t133 = -t164 * qJD(4) + t152 * t140;
t121 = t133 * pkin(5) + qJD(6) + t123;
t120 = -t133 * qJ(6) + t163;
t119 = t137 * pkin(5) - t134 * qJ(6) + t158;
t1 = [t155 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t149 ^ 2 / 0.2e1 + t151 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t155, -t142 * t160, t142 * t161 (-t135 * t148 + t136 * t150) * qJD(1), t136 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1 + t142 ^ 2 / 0.2e1, t140 ^ 2 / 0.2e1, -t140 * t139, t140 * qJD(4), -t139 * qJD(4), qJD(4) ^ 2 / 0.2e1, t157 * qJD(4) + t138 * t139, -t162 * qJD(4) + t138 * t140, t134 ^ 2 / 0.2e1, -t134 * t133, t134 * t137, -t133 * t137, t137 ^ 2 / 0.2e1, t123 * t133 + t158 * t137, t123 * t134 - t163 * t137, -t119 * t134 - t120 * t133, t120 ^ 2 / 0.2e1 + t119 ^ 2 / 0.2e1 + t121 ^ 2 / 0.2e1;];
T_reg  = t1;
