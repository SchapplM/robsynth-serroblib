% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% T_reg [1x23]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PPRRPR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:27
% EndTime: 2019-03-08 18:47:27
% DurationCPUTime: 0.16s
% Computational Cost: add. (282->42), mult. (703->96), div. (0->0), fcn. (580->14), ass. (0->45)
t198 = cos(pkin(6)) * qJD(1) + qJD(2);
t205 = sin(pkin(7));
t208 = cos(pkin(7));
t207 = cos(pkin(12));
t206 = sin(pkin(6));
t227 = qJD(1) * t206;
t221 = t207 * t227;
t233 = t198 * t205 + t208 * t221;
t211 = sin(qJ(3));
t214 = cos(qJ(3));
t204 = sin(pkin(12));
t222 = t204 * t227;
t232 = -t211 * t222 + t214 * t233;
t215 = qJD(3) ^ 2;
t231 = t215 / 0.2e1;
t230 = cos(pkin(13));
t223 = t233 * t211 + t214 * t222;
t186 = qJD(3) * pkin(9) + t223;
t190 = t208 * t198 - t205 * t221;
t210 = sin(qJ(4));
t213 = cos(qJ(4));
t228 = t213 * t186 + t210 * t190;
t179 = qJD(4) * qJ(5) + t228;
t182 = (-pkin(4) * t213 - qJ(5) * t210 - pkin(3)) * qJD(3) - t232;
t203 = sin(pkin(13));
t175 = t230 * t179 + t203 * t182;
t226 = qJD(3) * t210;
t225 = t213 * qJD(3);
t224 = qJD(3) * qJD(4);
t174 = -t203 * t179 + t230 * t182;
t220 = -t210 * t186 + t213 * t190;
t178 = -qJD(4) * pkin(4) + qJD(5) - t220;
t216 = qJD(1) ^ 2;
t212 = cos(qJ(6));
t209 = sin(qJ(6));
t199 = -qJD(6) + t225;
t194 = t203 * qJD(4) + t230 * t226;
t193 = -t230 * qJD(4) + t203 * t226;
t188 = -t209 * t193 + t212 * t194;
t187 = t212 * t193 + t209 * t194;
t185 = -qJD(3) * pkin(3) - t232;
t176 = t193 * pkin(5) + t178;
t173 = -t193 * pkin(10) + t175;
t172 = -pkin(5) * t225 - t194 * pkin(10) + t174;
t1 = [t216 / 0.2e1, t198 ^ 2 / 0.2e1 + (t204 ^ 2 / 0.2e1 + t207 ^ 2 / 0.2e1) * t216 * t206 ^ 2, t231, t232 * qJD(3), -t223 * qJD(3), t210 ^ 2 * t231, t213 * t215 * t210, t210 * t224, t213 * t224, qJD(4) ^ 2 / 0.2e1, t220 * qJD(4) - t185 * t225, -t228 * qJD(4) + t185 * t226, -t174 * t225 + t178 * t193, t175 * t225 + t178 * t194, -t174 * t194 - t175 * t193, t175 ^ 2 / 0.2e1 + t174 ^ 2 / 0.2e1 + t178 ^ 2 / 0.2e1, t188 ^ 2 / 0.2e1, -t188 * t187, -t188 * t199, t187 * t199, t199 ^ 2 / 0.2e1 -(t212 * t172 - t209 * t173) * t199 + t176 * t187 (t209 * t172 + t212 * t173) * t199 + t176 * t188;];
T_reg  = t1;
