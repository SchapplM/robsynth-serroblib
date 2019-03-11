% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRRRP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:01:15
% EndTime: 2019-03-09 06:01:15
% DurationCPUTime: 0.10s
% Computational Cost: add. (254->42), mult. (557->93), div. (0->0), fcn. (352->8), ass. (0->38)
t165 = qJD(1) ^ 2;
t178 = t165 / 0.2e1;
t177 = cos(qJ(5));
t161 = sin(qJ(4));
t163 = cos(qJ(4));
t162 = sin(qJ(3));
t173 = qJD(1) * t162;
t149 = t161 * qJD(3) + t163 * t173;
t164 = cos(qJ(3));
t172 = t164 * qJD(1);
t154 = -qJD(4) + t172;
t158 = sin(pkin(10));
t150 = (pkin(1) * t158 + pkin(7)) * qJD(1);
t174 = t162 * qJD(2) + t164 * t150;
t143 = qJD(3) * pkin(8) + t174;
t159 = cos(pkin(10));
t170 = -pkin(1) * t159 - pkin(2);
t144 = (-pkin(3) * t164 - pkin(8) * t162 + t170) * qJD(1);
t168 = -t161 * t143 + t163 * t144;
t132 = -t154 * pkin(4) - t149 * pkin(9) + t168;
t148 = -t163 * qJD(3) + t161 * t173;
t175 = t163 * t143 + t161 * t144;
t135 = -t148 * pkin(9) + t175;
t160 = sin(qJ(5));
t176 = t160 * t132 + t177 * t135;
t171 = qJD(1) * qJD(3);
t169 = t177 * t132 - t160 * t135;
t167 = t164 * qJD(2) - t162 * t150;
t142 = -qJD(3) * pkin(3) - t167;
t136 = t148 * pkin(4) + t142;
t152 = -qJD(5) + t154;
t151 = t170 * qJD(1);
t138 = -t160 * t148 + t177 * t149;
t137 = t177 * t148 + t160 * t149;
t133 = t137 * pkin(5) + qJD(6) + t136;
t129 = -t137 * qJ(6) + t176;
t128 = -t152 * pkin(5) - t138 * qJ(6) + t169;
t1 = [t178, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t158 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t165, t162 ^ 2 * t178, t164 * t165 * t162, t162 * t171, t164 * t171, qJD(3) ^ 2 / 0.2e1, t167 * qJD(3) - t151 * t172, -t174 * qJD(3) + t151 * t173, t149 ^ 2 / 0.2e1, -t149 * t148, -t149 * t154, t148 * t154, t154 ^ 2 / 0.2e1, t142 * t148 - t168 * t154, t142 * t149 + t175 * t154, t138 ^ 2 / 0.2e1, -t138 * t137, -t138 * t152, t137 * t152, t152 ^ 2 / 0.2e1, t136 * t137 - t169 * t152, t136 * t138 + t176 * t152, -t128 * t138 - t129 * t137, t129 ^ 2 / 0.2e1 + t128 ^ 2 / 0.2e1 + t133 ^ 2 / 0.2e1;];
T_reg  = t1;
