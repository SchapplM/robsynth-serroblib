% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:54:06
% EndTime: 2019-03-08 21:54:06
% DurationCPUTime: 0.13s
% Computational Cost: add. (312->47), mult. (768->102), div. (0->0), fcn. (596->12), ass. (0->44)
t204 = qJD(2) ^ 2;
t214 = t204 / 0.2e1;
t199 = sin(qJ(2));
t211 = qJD(1) * sin(pkin(6));
t186 = qJD(2) * pkin(8) + t199 * t211;
t202 = cos(qJ(3));
t210 = qJD(1) * cos(pkin(6));
t189 = t202 * t210;
t198 = sin(qJ(3));
t179 = qJD(3) * pkin(3) + t189 + (-qJ(4) * qJD(2) - t186) * t198;
t209 = qJD(2) * t202;
t212 = t202 * t186 + t198 * t210;
t180 = qJ(4) * t209 + t212;
t192 = sin(pkin(12));
t194 = cos(pkin(12));
t168 = t194 * t179 - t192 * t180;
t184 = (t192 * t202 + t194 * t198) * qJD(2);
t166 = qJD(3) * pkin(4) - t184 * pkin(9) + t168;
t169 = t192 * t179 + t194 * t180;
t183 = (-t192 * t198 + t194 * t202) * qJD(2);
t167 = t183 * pkin(9) + t169;
t197 = sin(qJ(5));
t201 = cos(qJ(5));
t213 = t197 * t166 + t201 * t167;
t208 = qJD(2) * qJD(3);
t203 = cos(qJ(2));
t207 = t203 * t211;
t206 = qJD(2) * t211;
t173 = -t201 * t183 + t197 * t184;
t205 = t201 * t166 - t197 * t167;
t182 = -t207 + qJD(4) + (-pkin(3) * t202 - pkin(2)) * qJD(2);
t175 = -t183 * pkin(4) + t182;
t200 = cos(qJ(6));
t196 = sin(qJ(6));
t191 = qJD(3) + qJD(5);
t187 = -qJD(2) * pkin(2) - t207;
t174 = t197 * t183 + t201 * t184;
t172 = qJD(6) + t173;
t171 = t200 * t174 + t196 * t191;
t170 = t196 * t174 - t200 * t191;
t163 = t173 * pkin(5) - t174 * pkin(10) + t175;
t162 = t191 * pkin(10) + t213;
t161 = -t191 * pkin(5) - t205;
t1 = [qJD(1) ^ 2 / 0.2e1, t214, t203 * t206, -t199 * t206, t198 ^ 2 * t214, t198 * t204 * t202, t198 * t208, t202 * t208, qJD(3) ^ 2 / 0.2e1 (-t198 * t186 + t189) * qJD(3) - t187 * t209, t187 * t198 * qJD(2) - t212 * qJD(3), -t168 * t184 + t169 * t183, t169 ^ 2 / 0.2e1 + t168 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1, t174 ^ 2 / 0.2e1, -t174 * t173, t174 * t191, -t173 * t191, t191 ^ 2 / 0.2e1, t175 * t173 + t205 * t191, t175 * t174 - t213 * t191, t171 ^ 2 / 0.2e1, -t171 * t170, t171 * t172, -t170 * t172, t172 ^ 2 / 0.2e1 (-t196 * t162 + t200 * t163) * t172 + t161 * t170 -(t200 * t162 + t196 * t163) * t172 + t161 * t171;];
T_reg  = t1;
