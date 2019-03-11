% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_energykin_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:15
% EndTime: 2019-03-09 09:15:15
% DurationCPUTime: 0.11s
% Computational Cost: add. (340->53), mult. (763->110), div. (0->0), fcn. (494->8), ass. (0->42)
t203 = qJD(1) ^ 2;
t213 = t203 / 0.2e1;
t202 = cos(qJ(2));
t212 = t202 * t203;
t199 = sin(qJ(2));
t210 = qJD(1) * t199;
t208 = pkin(7) * t210 + qJD(3);
t176 = -qJ(4) * t210 + (-pkin(2) - pkin(3)) * qJD(2) + t208;
t209 = qJD(1) * t202;
t184 = pkin(7) * t209 + qJD(2) * qJ(3);
t181 = -qJ(4) * t209 + t184;
t194 = sin(pkin(10));
t195 = cos(pkin(10));
t167 = t195 * t176 - t194 * t181;
t180 = (-t194 * t202 + t195 * t199) * qJD(1);
t163 = -qJD(2) * pkin(4) - t180 * pkin(8) + t167;
t168 = t194 * t176 + t195 * t181;
t179 = (t194 * t199 + t195 * t202) * qJD(1);
t164 = -t179 * pkin(8) + t168;
t198 = sin(qJ(5));
t201 = cos(qJ(5));
t211 = t198 * t163 + t201 * t164;
t207 = qJD(1) * qJD(2);
t182 = -qJD(1) * pkin(1) - pkin(2) * t209 - qJ(3) * t210;
t206 = t199 * t207;
t205 = t202 * t207;
t171 = t201 * t179 + t198 * t180;
t175 = pkin(3) * t209 + qJD(4) - t182;
t204 = t201 * t163 - t198 * t164;
t170 = t179 * pkin(4) + t175;
t200 = cos(qJ(6));
t197 = sin(qJ(6));
t191 = qJD(2) - qJD(5);
t183 = -qJD(2) * pkin(2) + t208;
t172 = -t198 * t179 + t201 * t180;
t169 = qJD(6) + t171;
t166 = t200 * t172 - t197 * t191;
t165 = t197 * t172 + t200 * t191;
t160 = t171 * pkin(5) - t172 * pkin(9) + t170;
t159 = -t191 * pkin(9) + t211;
t158 = t191 * pkin(5) - t204;
t1 = [t213, 0, 0, t199 ^ 2 * t213, t199 * t212, t206, t205, qJD(2) ^ 2 / 0.2e1, pkin(1) * t212 - pkin(7) * t206, -t203 * pkin(1) * t199 - pkin(7) * t205, -t183 * qJD(2) - t182 * t209 (t183 * t199 + t184 * t202) * qJD(1), t184 * qJD(2) - t182 * t210, t184 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1 + t183 ^ 2 / 0.2e1, -t167 * qJD(2) + t175 * t179, t168 * qJD(2) + t175 * t180, -t167 * t180 - t168 * t179, t168 ^ 2 / 0.2e1 + t167 ^ 2 / 0.2e1 + t175 ^ 2 / 0.2e1, t172 ^ 2 / 0.2e1, -t172 * t171, -t172 * t191, t171 * t191, t191 ^ 2 / 0.2e1, t170 * t171 - t204 * t191, t170 * t172 + t211 * t191, t166 ^ 2 / 0.2e1, -t166 * t165, t166 * t169, -t165 * t169, t169 ^ 2 / 0.2e1 (-t197 * t159 + t200 * t160) * t169 + t158 * t165 -(t200 * t159 + t197 * t160) * t169 + t158 * t166;];
T_reg  = t1;
