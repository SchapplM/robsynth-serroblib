% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:23:11
% EndTime: 2019-03-09 22:23:11
% DurationCPUTime: 0.15s
% Computational Cost: add. (650->53), mult. (1424->110), div. (0->0), fcn. (1067->10), ass. (0->48)
t206 = qJD(1) ^ 2;
t218 = t206 / 0.2e1;
t217 = cos(qJ(4));
t205 = cos(qJ(2));
t216 = t205 * t206;
t201 = sin(qJ(3));
t204 = cos(qJ(3));
t202 = sin(qJ(2));
t213 = qJD(1) * t202;
t184 = -t204 * qJD(2) + t201 * t213;
t185 = t201 * qJD(2) + t204 * t213;
t200 = sin(qJ(4));
t178 = -t200 * t184 + t217 * t185;
t212 = t205 * qJD(1);
t193 = -qJD(3) + t212;
t191 = -qJD(4) + t193;
t183 = (-pkin(2) * t205 - pkin(8) * t202 - pkin(1)) * qJD(1);
t190 = pkin(7) * t212 + qJD(2) * pkin(8);
t207 = t204 * t183 - t201 * t190;
t174 = -t193 * pkin(3) - t185 * pkin(9) + t207;
t214 = t201 * t183 + t204 * t190;
t176 = -t184 * pkin(9) + t214;
t208 = t217 * t174 - t200 * t176;
t164 = -t191 * pkin(4) - t178 * qJ(5) + t208;
t177 = t217 * t184 + t200 * t185;
t215 = t200 * t174 + t217 * t176;
t166 = -t177 * qJ(5) + t215;
t197 = sin(pkin(11));
t198 = cos(pkin(11));
t159 = t197 * t164 + t198 * t166;
t211 = qJD(1) * qJD(2);
t210 = t202 * t211;
t209 = t205 * t211;
t189 = -qJD(2) * pkin(2) + pkin(7) * t213;
t158 = t198 * t164 - t197 * t166;
t179 = t184 * pkin(3) + t189;
t171 = t177 * pkin(4) + qJD(5) + t179;
t203 = cos(qJ(6));
t199 = sin(qJ(6));
t187 = -qJD(6) + t191;
t170 = -t197 * t177 + t198 * t178;
t169 = -t198 * t177 - t197 * t178;
t167 = -t169 * pkin(5) + t171;
t163 = t199 * t169 + t203 * t170;
t162 = -t203 * t169 + t199 * t170;
t157 = t169 * pkin(10) + t159;
t156 = -t191 * pkin(5) - t170 * pkin(10) + t158;
t1 = [t218, 0, 0, t202 ^ 2 * t218, t202 * t216, t210, t209, qJD(2) ^ 2 / 0.2e1, pkin(1) * t216 - pkin(7) * t210, -t206 * pkin(1) * t202 - pkin(7) * t209, t185 ^ 2 / 0.2e1, -t185 * t184, -t185 * t193, t184 * t193, t193 ^ 2 / 0.2e1, t189 * t184 - t207 * t193, t189 * t185 + t214 * t193, t178 ^ 2 / 0.2e1, -t178 * t177, -t178 * t191, t177 * t191, t191 ^ 2 / 0.2e1, t179 * t177 - t208 * t191, t179 * t178 + t215 * t191, -t158 * t170 + t159 * t169, t159 ^ 2 / 0.2e1 + t158 ^ 2 / 0.2e1 + t171 ^ 2 / 0.2e1, t163 ^ 2 / 0.2e1, -t163 * t162, -t163 * t187, t162 * t187, t187 ^ 2 / 0.2e1 -(t203 * t156 - t199 * t157) * t187 + t167 * t162 (t199 * t156 + t203 * t157) * t187 + t167 * t163;];
T_reg  = t1;
