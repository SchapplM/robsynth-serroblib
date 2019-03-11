% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRR2
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
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:15
% EndTime: 2019-03-08 22:00:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (291->47), mult. (704->101), div. (0->0), fcn. (538->12), ass. (0->45)
t208 = qJD(2) ^ 2;
t220 = t208 / 0.2e1;
t219 = cos(qJ(5));
t204 = sin(qJ(2));
t216 = qJD(1) * sin(pkin(6));
t191 = qJD(2) * pkin(8) + t204 * t216;
t206 = cos(qJ(3));
t215 = qJD(1) * cos(pkin(6));
t195 = t206 * t215;
t203 = sin(qJ(3));
t180 = qJD(3) * pkin(3) + t195 + (-qJ(4) * qJD(2) - t191) * t203;
t213 = qJD(2) * t206;
t217 = t206 * t191 + t203 * t215;
t181 = qJ(4) * t213 + t217;
t197 = sin(pkin(12));
t199 = cos(pkin(12));
t171 = t197 * t180 + t199 * t181;
t169 = qJD(3) * pkin(9) + t171;
t207 = cos(qJ(2));
t211 = t207 * t216;
t186 = -t211 + qJD(4) + (-pkin(3) * t206 - pkin(2)) * qJD(2);
t214 = qJD(2) * t203;
t188 = -t197 * t214 + t199 * t213;
t189 = (t197 * t206 + t199 * t203) * qJD(2);
t176 = -t188 * pkin(4) - t189 * pkin(9) + t186;
t202 = sin(qJ(5));
t218 = t219 * t169 + t202 * t176;
t212 = qJD(2) * qJD(3);
t210 = qJD(2) * t216;
t209 = -t202 * t169 + t219 * t176;
t170 = t199 * t180 - t197 * t181;
t187 = qJD(5) - t188;
t168 = -qJD(3) * pkin(4) - t170;
t205 = cos(qJ(6));
t201 = sin(qJ(6));
t192 = -qJD(2) * pkin(2) - t211;
t185 = qJD(6) + t187;
t184 = t202 * qJD(3) + t219 * t189;
t183 = -t219 * qJD(3) + t202 * t189;
t173 = -t201 * t183 + t205 * t184;
t172 = t205 * t183 + t201 * t184;
t166 = t183 * pkin(5) + t168;
t165 = -t183 * pkin(10) + t218;
t164 = t187 * pkin(5) - t184 * pkin(10) + t209;
t1 = [qJD(1) ^ 2 / 0.2e1, t220, t207 * t210, -t204 * t210, t203 ^ 2 * t220, t203 * t208 * t206, t203 * t212, t206 * t212, qJD(3) ^ 2 / 0.2e1 (-t203 * t191 + t195) * qJD(3) - t192 * t213, -t217 * qJD(3) + t192 * t214, -t170 * t189 + t171 * t188, t171 ^ 2 / 0.2e1 + t170 ^ 2 / 0.2e1 + t186 ^ 2 / 0.2e1, t184 ^ 2 / 0.2e1, -t184 * t183, t184 * t187, -t183 * t187, t187 ^ 2 / 0.2e1, t168 * t183 + t209 * t187, t168 * t184 - t218 * t187, t173 ^ 2 / 0.2e1, -t173 * t172, t173 * t185, -t172 * t185, t185 ^ 2 / 0.2e1 (t205 * t164 - t201 * t165) * t185 + t166 * t172 -(t201 * t164 + t205 * t165) * t185 + t166 * t173;];
T_reg  = t1;
