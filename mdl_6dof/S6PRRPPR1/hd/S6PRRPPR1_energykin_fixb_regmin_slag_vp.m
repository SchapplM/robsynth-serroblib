% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:11
% EndTime: 2019-03-08 21:02:11
% DurationCPUTime: 0.13s
% Computational Cost: add. (355->48), mult. (843->102), div. (0->0), fcn. (630->12), ass. (0->44)
t209 = qJD(2) ^ 2;
t219 = t209 / 0.2e1;
t218 = cos(pkin(12));
t205 = sin(qJ(2));
t216 = qJD(1) * sin(pkin(6));
t192 = qJD(2) * pkin(8) + t205 * t216;
t207 = cos(qJ(3));
t215 = qJD(1) * cos(pkin(6));
t196 = t207 * t215;
t204 = sin(qJ(3));
t182 = qJD(3) * pkin(3) + t196 + (-qJ(4) * qJD(2) - t192) * t204;
t213 = qJD(2) * t207;
t217 = t207 * t192 + t204 * t215;
t183 = qJ(4) * t213 + t217;
t199 = sin(pkin(11));
t201 = cos(pkin(11));
t173 = t199 * t182 + t201 * t183;
t171 = qJD(3) * qJ(5) + t173;
t208 = cos(qJ(2));
t211 = t208 * t216;
t187 = -t211 + qJD(4) + (-pkin(3) * t207 - pkin(2)) * qJD(2);
t214 = qJD(2) * t204;
t189 = t199 * t214 - t201 * t213;
t190 = (t199 * t207 + t201 * t204) * qJD(2);
t178 = t189 * pkin(4) - t190 * qJ(5) + t187;
t198 = sin(pkin(12));
t167 = t218 * t171 + t198 * t178;
t212 = qJD(2) * qJD(3);
t210 = qJD(2) * t216;
t166 = -t198 * t171 + t218 * t178;
t172 = t201 * t182 - t199 * t183;
t170 = -qJD(3) * pkin(4) + qJD(5) - t172;
t206 = cos(qJ(6));
t203 = sin(qJ(6));
t193 = -qJD(2) * pkin(2) - t211;
t188 = qJD(6) + t189;
t186 = t198 * qJD(3) + t218 * t190;
t185 = -t218 * qJD(3) + t198 * t190;
t175 = -t203 * t185 + t206 * t186;
t174 = t206 * t185 + t203 * t186;
t168 = t185 * pkin(5) + t170;
t165 = -t185 * pkin(9) + t167;
t164 = t189 * pkin(5) - t186 * pkin(9) + t166;
t1 = [qJD(1) ^ 2 / 0.2e1, t219, t208 * t210, -t205 * t210, t204 ^ 2 * t219, t204 * t209 * t207, t204 * t212, t207 * t212, qJD(3) ^ 2 / 0.2e1 (-t204 * t192 + t196) * qJD(3) - t193 * t213, -t217 * qJD(3) + t193 * t214, -t172 * t190 - t173 * t189, t173 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1, t166 * t189 + t170 * t185, -t167 * t189 + t170 * t186, -t166 * t186 - t167 * t185, t167 ^ 2 / 0.2e1 + t166 ^ 2 / 0.2e1 + t170 ^ 2 / 0.2e1, t175 ^ 2 / 0.2e1, -t175 * t174, t175 * t188, -t174 * t188, t188 ^ 2 / 0.2e1 (t206 * t164 - t203 * t165) * t188 + t168 * t174 -(t203 * t164 + t206 * t165) * t188 + t168 * t175;];
T_reg  = t1;
