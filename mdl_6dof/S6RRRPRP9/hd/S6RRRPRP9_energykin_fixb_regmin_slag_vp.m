% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP9
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:27:14
% EndTime: 2019-03-09 17:27:14
% DurationCPUTime: 0.12s
% Computational Cost: add. (385->49), mult. (768->99), div. (0->0), fcn. (485->6), ass. (0->38)
t184 = qJD(1) ^ 2;
t196 = t184 / 0.2e1;
t183 = cos(qJ(2));
t195 = t183 * t184;
t179 = sin(qJ(3));
t182 = cos(qJ(3));
t180 = sin(qJ(2));
t192 = qJD(1) * t180;
t166 = t179 * qJD(2) + t182 * t192;
t191 = qJD(1) * t183;
t173 = -qJD(3) + t191;
t162 = (-pkin(2) * t183 - pkin(8) * t180 - pkin(1)) * qJD(1);
t170 = pkin(7) * t191 + qJD(2) * pkin(8);
t187 = t162 * t182 - t179 * t170;
t186 = qJD(4) - t187;
t151 = -t166 * pkin(9) + (pkin(3) + pkin(4)) * t173 + t186;
t193 = t179 * t162 + t182 * t170;
t156 = -t173 * qJ(4) + t193;
t165 = -t182 * qJD(2) + t179 * t192;
t153 = pkin(9) * t165 + t156;
t178 = sin(qJ(5));
t181 = cos(qJ(5));
t194 = t178 * t151 + t181 * t153;
t169 = -qJD(2) * pkin(2) + pkin(7) * t192;
t190 = qJD(1) * qJD(2);
t189 = t180 * t190;
t188 = t183 * t190;
t157 = t165 * pkin(3) - t166 * qJ(4) + t169;
t185 = t151 * t181 - t178 * t153;
t154 = -pkin(4) * t165 - t157;
t172 = qJD(5) + t173;
t159 = t165 * t178 + t166 * t181;
t158 = -t181 * t165 + t166 * t178;
t155 = t173 * pkin(3) + t186;
t149 = pkin(5) * t158 - qJ(6) * t159 + t154;
t148 = qJ(6) * t172 + t194;
t147 = -pkin(5) * t172 + qJD(6) - t185;
t1 = [t196, 0, 0, t180 ^ 2 * t196, t180 * t195, t189, t188, qJD(2) ^ 2 / 0.2e1, pkin(1) * t195 - pkin(7) * t189, -pkin(1) * t180 * t184 - pkin(7) * t188, t166 ^ 2 / 0.2e1, -t166 * t165, -t166 * t173, t165 * t173, t173 ^ 2 / 0.2e1, t169 * t165 - t187 * t173, t169 * t166 + t193 * t173, t155 * t173 + t157 * t165, t155 * t166 - t156 * t165, -t156 * t173 - t157 * t166, t156 ^ 2 / 0.2e1 + t157 ^ 2 / 0.2e1 + t155 ^ 2 / 0.2e1, t159 ^ 2 / 0.2e1, -t159 * t158, t159 * t172, -t158 * t172, t172 ^ 2 / 0.2e1, t154 * t158 + t185 * t172, t154 * t159 - t194 * t172, -t147 * t172 + t149 * t158, t147 * t159 - t148 * t158, t148 * t172 - t149 * t159, t148 ^ 2 / 0.2e1 + t149 ^ 2 / 0.2e1 + t147 ^ 2 / 0.2e1;];
T_reg  = t1;
