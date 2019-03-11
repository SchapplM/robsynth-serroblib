% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:45:55
% EndTime: 2019-03-09 16:45:55
% DurationCPUTime: 0.13s
% Computational Cost: add. (382->48), mult. (811->98), div. (0->0), fcn. (534->6), ass. (0->41)
t183 = pkin(3) + pkin(9);
t182 = -pkin(8) - pkin(7);
t168 = qJD(1) ^ 2;
t181 = t168 / 0.2e1;
t167 = cos(qJ(2));
t180 = t167 * t168;
t163 = sin(qJ(3));
t164 = sin(qJ(2));
t166 = cos(qJ(3));
t153 = (t163 * t167 + t164 * t166) * qJD(1);
t161 = qJD(2) + qJD(3);
t177 = qJD(1) * t164;
t156 = qJD(2) * pkin(2) + t182 * t177;
t176 = qJD(1) * t167;
t157 = t182 * t176;
t172 = t166 * t156 + t163 * t157;
t171 = qJD(4) - t172;
t141 = t153 * pkin(4) - t183 * t161 + t171;
t152 = t163 * t177 - t166 * t176;
t158 = (-pkin(2) * t167 - pkin(1)) * qJD(1);
t169 = -t153 * qJ(4) + t158;
t142 = t183 * t152 + t169;
t162 = sin(qJ(5));
t165 = cos(qJ(5));
t179 = t162 * t141 + t165 * t142;
t178 = t163 * t156 - t166 * t157;
t175 = qJD(1) * qJD(2);
t146 = -t161 * qJ(4) - t178;
t174 = t164 * t175;
t173 = t167 * t175;
t170 = t165 * t141 - t162 * t142;
t143 = -t152 * pkin(4) - t146;
t151 = qJD(5) + t153;
t148 = t162 * t152 + t165 * t161;
t147 = -t165 * t152 + t162 * t161;
t145 = -t161 * pkin(3) + t171;
t144 = t152 * pkin(3) + t169;
t138 = t147 * pkin(5) - t148 * qJ(6) + t143;
t137 = t151 * qJ(6) + t179;
t136 = -t151 * pkin(5) + qJD(6) - t170;
t1 = [t181, 0, 0, t164 ^ 2 * t181, t164 * t180, t174, t173, qJD(2) ^ 2 / 0.2e1, pkin(1) * t180 - pkin(7) * t174, -t168 * pkin(1) * t164 - pkin(7) * t173, t153 ^ 2 / 0.2e1, -t153 * t152, t153 * t161, -t152 * t161, t161 ^ 2 / 0.2e1, t158 * t152 + t172 * t161, t158 * t153 - t178 * t161, t145 * t153 + t146 * t152, -t144 * t152 + t145 * t161, -t144 * t153 - t146 * t161, t144 ^ 2 / 0.2e1 + t146 ^ 2 / 0.2e1 + t145 ^ 2 / 0.2e1, t148 ^ 2 / 0.2e1, -t148 * t147, t148 * t151, -t147 * t151, t151 ^ 2 / 0.2e1, t143 * t147 + t170 * t151, t143 * t148 - t179 * t151, -t136 * t151 + t138 * t147, t136 * t148 - t137 * t147, t137 * t151 - t138 * t148, t137 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1;];
T_reg  = t1;
