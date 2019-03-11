% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:51:39
% EndTime: 2019-03-09 20:51:39
% DurationCPUTime: 0.15s
% Computational Cost: add. (455->48), mult. (951->98), div. (0->0), fcn. (649->6), ass. (0->41)
t185 = -pkin(4) - pkin(5);
t184 = -pkin(8) - pkin(7);
t168 = qJD(1) ^ 2;
t183 = t168 / 0.2e1;
t182 = cos(qJ(4));
t167 = cos(qJ(2));
t181 = t167 * t168;
t164 = sin(qJ(3));
t166 = cos(qJ(3));
t176 = qJD(1) * t167;
t165 = sin(qJ(2));
t177 = qJD(1) * t165;
t153 = t164 * t177 - t166 * t176;
t154 = (t164 * t167 + t165 * t166) * qJD(1);
t159 = (-pkin(2) * t167 - pkin(1)) * qJD(1);
t142 = t153 * pkin(3) - t154 * pkin(9) + t159;
t162 = qJD(2) + qJD(3);
t157 = qJD(2) * pkin(2) + t184 * t177;
t158 = t184 * t176;
t178 = t164 * t157 - t166 * t158;
t146 = t162 * pkin(9) + t178;
t163 = sin(qJ(4));
t180 = t163 * t142 + t182 * t146;
t179 = t166 * t157 + t164 * t158;
t175 = qJD(1) * qJD(2);
t152 = qJD(4) + t153;
t139 = t152 * qJ(5) + t180;
t174 = t165 * t175;
t173 = t167 * t175;
t172 = t162 * pkin(3) + t179;
t171 = t182 * t142 - t163 * t146;
t170 = qJD(5) - t171;
t148 = t182 * t154 + t163 * t162;
t169 = t148 * qJ(5) + t172;
t147 = t163 * t154 - t182 * t162;
t140 = t147 * pkin(4) - t169;
t138 = -t152 * pkin(4) + t170;
t137 = t185 * t147 + qJD(6) + t169;
t136 = t147 * qJ(6) + t139;
t135 = -t148 * qJ(6) + t185 * t152 + t170;
t1 = [t183, 0, 0, t165 ^ 2 * t183, t165 * t181, t174, t173, qJD(2) ^ 2 / 0.2e1, pkin(1) * t181 - pkin(7) * t174, -t168 * pkin(1) * t165 - pkin(7) * t173, t154 ^ 2 / 0.2e1, -t154 * t153, t154 * t162, -t153 * t162, t162 ^ 2 / 0.2e1, t159 * t153 + t179 * t162, t159 * t154 - t178 * t162, t148 ^ 2 / 0.2e1, -t148 * t147, t148 * t152, -t147 * t152, t152 ^ 2 / 0.2e1, -t147 * t172 + t171 * t152, -t148 * t172 - t180 * t152, -t138 * t152 + t140 * t147, t138 * t148 - t139 * t147, t139 * t152 - t140 * t148, t139 ^ 2 / 0.2e1 + t140 ^ 2 / 0.2e1 + t138 ^ 2 / 0.2e1, -t135 * t152 - t137 * t147, t136 * t152 + t137 * t148, -t135 * t148 + t136 * t147, t136 ^ 2 / 0.2e1 + t135 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1;];
T_reg  = t1;
