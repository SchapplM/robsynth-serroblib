% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRR5
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
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:20:46
% EndTime: 2019-03-08 22:20:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (335->48), mult. (774->103), div. (0->0), fcn. (584->12), ass. (0->45)
t206 = qJD(2) ^ 2;
t219 = t206 / 0.2e1;
t218 = cos(qJ(5));
t202 = sin(qJ(2));
t215 = qJD(1) * sin(pkin(6));
t187 = qJD(2) * pkin(8) + t202 * t215;
t201 = sin(qJ(3));
t204 = cos(qJ(3));
t214 = qJD(1) * cos(pkin(6));
t216 = t204 * t187 + t201 * t214;
t180 = qJD(3) * qJ(4) + t216;
t205 = cos(qJ(2));
t210 = t205 * t215;
t181 = -t210 + (-pkin(3) * t204 - qJ(4) * t201 - pkin(2)) * qJD(2);
t195 = sin(pkin(12));
t197 = cos(pkin(12));
t170 = -t195 * t180 + t197 * t181;
t213 = qJD(2) * t201;
t186 = t195 * qJD(3) + t197 * t213;
t212 = t204 * qJD(2);
t165 = -pkin(4) * t212 - t186 * pkin(9) + t170;
t171 = t197 * t180 + t195 * t181;
t185 = -t197 * qJD(3) + t195 * t213;
t169 = -t185 * pkin(9) + t171;
t200 = sin(qJ(5));
t217 = t200 * t165 + t218 * t169;
t211 = qJD(2) * qJD(3);
t209 = qJD(2) * t215;
t208 = t218 * t165 - t200 * t169;
t192 = -qJD(5) + t212;
t207 = -t201 * t187 + t204 * t214;
t177 = -qJD(3) * pkin(3) + qJD(4) - t207;
t172 = t185 * pkin(4) + t177;
t203 = cos(qJ(6));
t199 = sin(qJ(6));
t189 = -qJD(6) + t192;
t188 = -qJD(2) * pkin(2) - t210;
t175 = -t200 * t185 + t218 * t186;
t174 = t218 * t185 + t200 * t186;
t167 = -t199 * t174 + t203 * t175;
t166 = t203 * t174 + t199 * t175;
t164 = t174 * pkin(5) + t172;
t161 = -t174 * pkin(10) + t217;
t160 = -t192 * pkin(5) - t175 * pkin(10) + t208;
t1 = [qJD(1) ^ 2 / 0.2e1, t219, t205 * t209, -t202 * t209, t201 ^ 2 * t219, t201 * t206 * t204, t201 * t211, t204 * t211, qJD(3) ^ 2 / 0.2e1, t207 * qJD(3) - t188 * t212, -t216 * qJD(3) + t188 * t213, -t170 * t212 + t177 * t185, t171 * t212 + t177 * t186, -t170 * t186 - t171 * t185, t171 ^ 2 / 0.2e1 + t170 ^ 2 / 0.2e1 + t177 ^ 2 / 0.2e1, t175 ^ 2 / 0.2e1, -t175 * t174, -t175 * t192, t174 * t192, t192 ^ 2 / 0.2e1, t172 * t174 - t208 * t192, t172 * t175 + t217 * t192, t167 ^ 2 / 0.2e1, -t167 * t166, -t167 * t189, t166 * t189, t189 ^ 2 / 0.2e1 -(t203 * t160 - t199 * t161) * t189 + t164 * t166 (t199 * t160 + t203 * t161) * t189 + t164 * t167;];
T_reg  = t1;
