% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPP3
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:56:28
% EndTime: 2019-03-09 20:56:28
% DurationCPUTime: 0.14s
% Computational Cost: add. (455->48), mult. (951->98), div. (0->0), fcn. (649->6), ass. (0->41)
t180 = -pkin(8) - pkin(7);
t165 = qJD(1) ^ 2;
t179 = t165 / 0.2e1;
t178 = pkin(4) + qJ(6);
t164 = cos(qJ(2));
t177 = t164 * t165;
t160 = sin(qJ(3));
t163 = cos(qJ(3));
t173 = qJD(1) * t164;
t161 = sin(qJ(2));
t174 = qJD(1) * t161;
t149 = t160 * t174 - t163 * t173;
t150 = (t160 * t164 + t161 * t163) * qJD(1);
t155 = (-pkin(2) * t164 - pkin(1)) * qJD(1);
t140 = t149 * pkin(3) - t150 * pkin(9) + t155;
t158 = qJD(2) + qJD(3);
t153 = qJD(2) * pkin(2) + t180 * t174;
t154 = t180 * t173;
t175 = t160 * t153 - t163 * t154;
t144 = t158 * pkin(9) + t175;
t159 = sin(qJ(4));
t162 = cos(qJ(4));
t176 = t159 * t140 + t162 * t144;
t172 = qJD(1) * qJD(2);
t171 = t161 * t172;
t170 = t164 * t172;
t169 = t162 * t140 - t159 * t144;
t168 = t163 * t153 + t160 * t154;
t148 = qJD(4) + t149;
t137 = -t148 * qJ(5) - t176;
t167 = qJD(5) - t169;
t143 = -t158 * pkin(3) - t168;
t146 = t162 * t150 + t159 * t158;
t166 = -t146 * qJ(5) + t143;
t145 = t159 * t150 - t162 * t158;
t138 = t145 * pkin(4) + t166;
t136 = -t148 * pkin(4) + t167;
t135 = t178 * t145 + t166;
t134 = -t145 * pkin(5) + qJD(6) - t137;
t133 = t146 * pkin(5) - t178 * t148 + t167;
t1 = [t179, 0, 0, t161 ^ 2 * t179, t161 * t177, t171, t170, qJD(2) ^ 2 / 0.2e1, pkin(1) * t177 - pkin(7) * t171, -t165 * pkin(1) * t161 - pkin(7) * t170, t150 ^ 2 / 0.2e1, -t150 * t149, t150 * t158, -t149 * t158, t158 ^ 2 / 0.2e1, t155 * t149 + t168 * t158, t155 * t150 - t175 * t158, t146 ^ 2 / 0.2e1, -t146 * t145, t146 * t148, -t145 * t148, t148 ^ 2 / 0.2e1, t143 * t145 + t169 * t148, t143 * t146 - t176 * t148, t136 * t146 + t137 * t145, t136 * t148 - t138 * t145, -t137 * t148 - t138 * t146, t138 ^ 2 / 0.2e1 + t137 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1, t133 * t146 - t134 * t145, t134 * t148 - t135 * t146, -t133 * t148 + t135 * t145, t135 ^ 2 / 0.2e1 + t133 ^ 2 / 0.2e1 + t134 ^ 2 / 0.2e1;];
T_reg  = t1;
