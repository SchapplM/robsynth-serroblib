% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:36
% EndTime: 2019-03-09 01:42:36
% DurationCPUTime: 0.12s
% Computational Cost: add. (220->48), mult. (522->96), div. (0->0), fcn. (329->8), ass. (0->37)
t166 = pkin(4) + pkin(8);
t150 = sin(pkin(9));
t145 = (pkin(1) * t150 + qJ(3)) * qJD(1);
t151 = cos(pkin(10));
t148 = t151 * qJD(2);
t149 = sin(pkin(10));
t132 = t148 + (-pkin(7) * qJD(1) - t145) * t149;
t137 = t149 * qJD(2) + t151 * t145;
t163 = qJD(1) * t151;
t133 = pkin(7) * t163 + t137;
t154 = sin(qJ(4));
t156 = cos(qJ(4));
t165 = t154 * t132 + t156 * t133;
t164 = qJD(1) * t149;
t152 = cos(pkin(9));
t162 = -pkin(1) * t152 - pkin(2);
t161 = t132 * t156 - t154 * t133;
t127 = -qJD(4) * qJ(5) - t165;
t160 = qJD(5) - t161;
t142 = (t149 * t156 + t151 * t154) * qJD(1);
t140 = qJD(3) + (-pkin(3) * t151 + t162) * qJD(1);
t159 = -qJ(5) * t142 + t140;
t157 = qJD(1) ^ 2;
t155 = cos(qJ(6));
t153 = sin(qJ(6));
t144 = t162 * qJD(1) + qJD(3);
t141 = t154 * t164 - t156 * t163;
t139 = qJD(6) + t142;
t136 = -t145 * t149 + t148;
t135 = qJD(4) * t155 + t141 * t153;
t134 = qJD(4) * t153 - t155 * t141;
t128 = pkin(4) * t141 + t159;
t126 = -qJD(4) * pkin(4) + t160;
t125 = t166 * t141 + t159;
t124 = -pkin(5) * t141 - t127;
t123 = pkin(5) * t142 - t166 * qJD(4) + t160;
t1 = [t157 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t150 ^ 2 / 0.2e1 + t152 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t157, -t144 * t163, t144 * t164 (-t136 * t149 + t137 * t151) * qJD(1), t137 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1 + t144 ^ 2 / 0.2e1, t142 ^ 2 / 0.2e1, -t142 * t141, t142 * qJD(4), -t141 * qJD(4), qJD(4) ^ 2 / 0.2e1, t161 * qJD(4) + t140 * t141, -t165 * qJD(4) + t140 * t142, t126 * t142 + t127 * t141, qJD(4) * t126 - t128 * t141, -qJD(4) * t127 - t128 * t142, t128 ^ 2 / 0.2e1 + t127 ^ 2 / 0.2e1 + t126 ^ 2 / 0.2e1, t135 ^ 2 / 0.2e1, -t135 * t134, t135 * t139, -t134 * t139, t139 ^ 2 / 0.2e1 (t123 * t155 - t153 * t125) * t139 + t124 * t134 -(t123 * t153 + t125 * t155) * t139 + t124 * t135;];
T_reg  = t1;
