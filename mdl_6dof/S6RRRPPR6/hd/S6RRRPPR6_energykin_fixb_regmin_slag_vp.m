% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:53:54
% EndTime: 2019-03-09 15:53:54
% DurationCPUTime: 0.13s
% Computational Cost: add. (632->54), mult. (1511->109), div. (0->0), fcn. (1165->10), ass. (0->48)
t236 = pkin(4) + pkin(10);
t235 = cos(qJ(3));
t213 = sin(pkin(6));
t221 = qJD(1) ^ 2;
t234 = t213 ^ 2 * t221;
t230 = cos(pkin(6)) * qJD(1);
t210 = qJD(2) + t230;
t217 = sin(qJ(3));
t218 = sin(qJ(2));
t231 = qJD(1) * t213;
t227 = t218 * t231;
t203 = t217 * t210 + t235 * t227;
t220 = cos(qJ(2));
t226 = t220 * t231;
t205 = -qJD(3) + t226;
t229 = pkin(1) * t230;
t232 = pkin(8) * t226 + t218 * t229;
t199 = t210 * pkin(9) + t232;
t201 = (-pkin(2) * t220 - pkin(9) * t218 - pkin(1)) * t231;
t225 = -t217 * t199 + t235 * t201;
t184 = -t205 * pkin(3) - t203 * qJ(4) + t225;
t202 = -t235 * t210 + t217 * t227;
t233 = t235 * t199 + t217 * t201;
t187 = -t202 * qJ(4) + t233;
t212 = sin(pkin(11));
t214 = cos(pkin(11));
t181 = t212 * t184 + t214 * t187;
t228 = t220 * t234;
t180 = t214 * t184 - t212 * t187;
t179 = t205 * qJ(5) - t181;
t224 = qJD(5) - t180;
t194 = -t212 * t202 + t214 * t203;
t223 = -pkin(8) * t227 + t220 * t229;
t198 = -t210 * pkin(2) - t223;
t192 = t202 * pkin(3) + qJD(4) + t198;
t222 = -t194 * qJ(5) + t192;
t219 = cos(qJ(6));
t216 = sin(qJ(6));
t193 = t214 * t202 + t212 * t203;
t191 = qJD(6) + t194;
t189 = t216 * t193 - t219 * t205;
t188 = -t219 * t193 - t216 * t205;
t182 = t193 * pkin(4) + t222;
t178 = t205 * pkin(4) + t224;
t177 = t236 * t193 + t222;
t176 = -t193 * pkin(5) - t179;
t175 = t194 * pkin(5) + t236 * t205 + t224;
t1 = [t221 / 0.2e1, 0, 0, t218 ^ 2 * t234 / 0.2e1, t218 * t228, t210 * t227, t210 * t226, t210 ^ 2 / 0.2e1, pkin(1) * t228 + t223 * t210, -pkin(1) * t218 * t234 - t232 * t210, t203 ^ 2 / 0.2e1, -t203 * t202, -t203 * t205, t202 * t205, t205 ^ 2 / 0.2e1, t198 * t202 - t225 * t205, t198 * t203 + t233 * t205, -t180 * t194 - t181 * t193, t181 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1 + t192 ^ 2 / 0.2e1, t178 * t194 + t179 * t193, -t178 * t205 - t182 * t193, t179 * t205 - t182 * t194, t182 ^ 2 / 0.2e1 + t179 ^ 2 / 0.2e1 + t178 ^ 2 / 0.2e1, t189 ^ 2 / 0.2e1, -t189 * t188, t189 * t191, -t188 * t191, t191 ^ 2 / 0.2e1 (t219 * t175 - t216 * t177) * t191 + t176 * t188 -(t216 * t175 + t219 * t177) * t191 + t176 * t189;];
T_reg  = t1;
