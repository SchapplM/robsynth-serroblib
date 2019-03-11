% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:35:47
% EndTime: 2019-03-10 03:35:47
% DurationCPUTime: 0.16s
% Computational Cost: add. (583->53), mult. (1305->112), div. (0->0), fcn. (1019->10), ass. (0->50)
t212 = -pkin(8) - pkin(7);
t196 = qJD(1) ^ 2;
t211 = t196 / 0.2e1;
t210 = cos(qJ(3));
t209 = cos(qJ(5));
t195 = cos(qJ(2));
t208 = t195 * t196;
t187 = qJD(2) + qJD(3);
t186 = qJD(4) + t187;
t191 = sin(qJ(3));
t192 = sin(qJ(2));
t179 = (t191 * t195 + t210 * t192) * qJD(1);
t204 = qJD(1) * t192;
t181 = qJD(2) * pkin(2) + t212 * t204;
t203 = qJD(1) * t195;
t182 = t212 * t203;
t197 = t210 * t181 + t191 * t182;
t163 = t187 * pkin(3) - t179 * pkin(9) + t197;
t178 = t191 * t204 - t210 * t203;
t205 = t191 * t181 - t210 * t182;
t167 = -t178 * pkin(9) + t205;
t190 = sin(qJ(4));
t194 = cos(qJ(4));
t206 = t190 * t163 + t194 * t167;
t156 = t186 * pkin(10) + t206;
t172 = t194 * t178 + t190 * t179;
t173 = -t190 * t178 + t194 * t179;
t183 = (-pkin(2) * t195 - pkin(1)) * qJD(1);
t174 = t178 * pkin(3) + t183;
t161 = t172 * pkin(4) - t173 * pkin(10) + t174;
t189 = sin(qJ(5));
t207 = t209 * t156 + t189 * t161;
t202 = qJD(1) * qJD(2);
t201 = t192 * t202;
t200 = t195 * t202;
t199 = -t189 * t156 + t209 * t161;
t198 = t194 * t163 - t190 * t167;
t171 = qJD(5) + t172;
t155 = -t186 * pkin(4) - t198;
t193 = cos(qJ(6));
t188 = sin(qJ(6));
t170 = qJD(6) + t171;
t169 = t209 * t173 + t189 * t186;
t168 = t189 * t173 - t209 * t186;
t160 = -t188 * t168 + t193 * t169;
t159 = t193 * t168 + t188 * t169;
t153 = t168 * pkin(5) + t155;
t152 = -t168 * pkin(11) + t207;
t151 = t171 * pkin(5) - t169 * pkin(11) + t199;
t1 = [t211, 0, 0, t192 ^ 2 * t211, t192 * t208, t201, t200, qJD(2) ^ 2 / 0.2e1, pkin(1) * t208 - pkin(7) * t201, -t196 * pkin(1) * t192 - pkin(7) * t200, t179 ^ 2 / 0.2e1, -t179 * t178, t179 * t187, -t178 * t187, t187 ^ 2 / 0.2e1, t183 * t178 + t197 * t187, t183 * t179 - t205 * t187, t173 ^ 2 / 0.2e1, -t173 * t172, t173 * t186, -t172 * t186, t186 ^ 2 / 0.2e1, t174 * t172 + t198 * t186, t174 * t173 - t206 * t186, t169 ^ 2 / 0.2e1, -t169 * t168, t169 * t171, -t168 * t171, t171 ^ 2 / 0.2e1, t155 * t168 + t199 * t171, t155 * t169 - t207 * t171, t160 ^ 2 / 0.2e1, -t160 * t159, t160 * t170, -t159 * t170, t170 ^ 2 / 0.2e1 (t193 * t151 - t188 * t152) * t170 + t153 * t159 -(t188 * t151 + t193 * t152) * t170 + t153 * t160;];
T_reg  = t1;
