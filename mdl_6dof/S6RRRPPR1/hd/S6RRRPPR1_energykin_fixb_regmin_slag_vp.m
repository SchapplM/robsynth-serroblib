% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:23:05
% EndTime: 2019-03-09 15:23:05
% DurationCPUTime: 0.16s
% Computational Cost: add. (665->53), mult. (1557->110), div. (0->0), fcn. (1160->10), ass. (0->48)
t215 = -pkin(8) - pkin(7);
t203 = qJD(1) ^ 2;
t214 = t203 / 0.2e1;
t213 = cos(qJ(3));
t212 = cos(pkin(11));
t202 = cos(qJ(2));
t211 = t202 * t203;
t199 = sin(qJ(3));
t200 = sin(qJ(2));
t187 = (t199 * t202 + t213 * t200) * qJD(1);
t194 = qJD(2) + qJD(3);
t209 = qJD(1) * t200;
t189 = qJD(2) * pkin(2) + t215 * t209;
t208 = qJD(1) * t202;
t190 = t215 * t208;
t204 = t213 * t189 + t199 * t190;
t172 = t194 * pkin(3) - t187 * qJ(4) + t204;
t186 = t199 * t209 - t213 * t208;
t210 = t199 * t189 - t213 * t190;
t176 = -t186 * qJ(4) + t210;
t196 = sin(pkin(10));
t197 = cos(pkin(10));
t165 = t196 * t172 + t197 * t176;
t163 = t194 * qJ(5) + t165;
t180 = t197 * t186 + t196 * t187;
t181 = -t196 * t186 + t197 * t187;
t191 = (-pkin(2) * t202 - pkin(1)) * qJD(1);
t182 = t186 * pkin(3) + qJD(4) + t191;
t168 = t180 * pkin(4) - t181 * qJ(5) + t182;
t195 = sin(pkin(11));
t159 = t212 * t163 + t195 * t168;
t207 = qJD(1) * qJD(2);
t206 = t200 * t207;
t205 = t202 * t207;
t158 = -t195 * t163 + t212 * t168;
t164 = t197 * t172 - t196 * t176;
t162 = -t194 * pkin(4) + qJD(5) - t164;
t201 = cos(qJ(6));
t198 = sin(qJ(6));
t179 = qJD(6) + t180;
t178 = t212 * t181 + t195 * t194;
t177 = t195 * t181 - t212 * t194;
t170 = -t198 * t177 + t201 * t178;
t169 = t201 * t177 + t198 * t178;
t160 = t177 * pkin(5) + t162;
t157 = -t177 * pkin(9) + t159;
t156 = t180 * pkin(5) - t178 * pkin(9) + t158;
t1 = [t214, 0, 0, t200 ^ 2 * t214, t200 * t211, t206, t205, qJD(2) ^ 2 / 0.2e1, pkin(1) * t211 - pkin(7) * t206, -t203 * pkin(1) * t200 - pkin(7) * t205, t187 ^ 2 / 0.2e1, -t187 * t186, t187 * t194, -t186 * t194, t194 ^ 2 / 0.2e1, t191 * t186 + t204 * t194, t191 * t187 - t210 * t194, -t164 * t181 - t165 * t180, t165 ^ 2 / 0.2e1 + t164 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1, t158 * t180 + t162 * t177, -t159 * t180 + t162 * t178, -t158 * t178 - t159 * t177, t159 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1 + t162 ^ 2 / 0.2e1, t170 ^ 2 / 0.2e1, -t170 * t169, t170 * t179, -t169 * t179, t179 ^ 2 / 0.2e1 (t201 * t156 - t198 * t157) * t179 + t160 * t169 -(t198 * t156 + t201 * t157) * t179 + t160 * t170;];
T_reg  = t1;
