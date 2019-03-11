% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:35:59
% EndTime: 2019-03-09 15:35:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (444->51), mult. (964->103), div. (0->0), fcn. (662->8), ass. (0->44)
t211 = pkin(4) + pkin(5);
t197 = qJD(1) ^ 2;
t210 = t197 / 0.2e1;
t209 = cos(qJ(3));
t196 = cos(qJ(2));
t208 = t196 * t197;
t193 = sin(qJ(3));
t194 = sin(qJ(2));
t206 = qJD(1) * t194;
t180 = t193 * qJD(2) + t209 * t206;
t205 = t196 * qJD(1);
t186 = -qJD(3) + t205;
t178 = (-pkin(2) * t196 - pkin(8) * t194 - pkin(1)) * qJD(1);
t183 = pkin(7) * t205 + qJD(2) * pkin(8);
t201 = t209 * t178 - t193 * t183;
t167 = -t186 * pkin(3) - t180 * qJ(4) + t201;
t179 = -t209 * qJD(2) + t193 * t206;
t207 = t193 * t178 + t209 * t183;
t170 = -t179 * qJ(4) + t207;
t190 = sin(pkin(10));
t191 = cos(pkin(10));
t162 = t190 * t167 + t191 * t170;
t204 = qJD(1) * qJD(2);
t160 = -t186 * qJ(5) + t162;
t203 = t194 * t204;
t202 = t196 * t204;
t182 = -qJD(2) * pkin(2) + pkin(7) * t206;
t161 = t191 * t167 - t190 * t170;
t200 = qJD(5) - t161;
t199 = -t179 * pkin(3) - qJD(4) - t182;
t173 = -t190 * t179 + t191 * t180;
t198 = t173 * qJ(5) + t199;
t195 = cos(qJ(6));
t192 = sin(qJ(6));
t185 = qJD(6) + t186;
t172 = t191 * t179 + t190 * t180;
t165 = t192 * t172 + t195 * t173;
t164 = -t195 * t172 + t192 * t173;
t163 = t172 * pkin(4) - t198;
t159 = t186 * pkin(4) + t200;
t158 = -t211 * t172 + t198;
t157 = t172 * pkin(9) + t160;
t156 = -t173 * pkin(9) + t211 * t186 + t200;
t1 = [t210, 0, 0, t194 ^ 2 * t210, t194 * t208, t203, t202, qJD(2) ^ 2 / 0.2e1, pkin(1) * t208 - pkin(7) * t203, -t197 * pkin(1) * t194 - pkin(7) * t202, t180 ^ 2 / 0.2e1, -t180 * t179, -t180 * t186, t179 * t186, t186 ^ 2 / 0.2e1, t182 * t179 - t201 * t186, t182 * t180 + t207 * t186, -t161 * t173 - t162 * t172, t162 ^ 2 / 0.2e1 + t161 ^ 2 / 0.2e1 + t199 ^ 2 / 0.2e1, t159 * t186 + t163 * t172, t159 * t173 - t160 * t172, -t160 * t186 - t163 * t173, t160 ^ 2 / 0.2e1 + t163 ^ 2 / 0.2e1 + t159 ^ 2 / 0.2e1, t165 ^ 2 / 0.2e1, -t165 * t164, t165 * t185, -t164 * t185, t185 ^ 2 / 0.2e1 (t195 * t156 - t192 * t157) * t185 + t158 * t164 -(t192 * t156 + t195 * t157) * t185 + t158 * t165;];
T_reg  = t1;
