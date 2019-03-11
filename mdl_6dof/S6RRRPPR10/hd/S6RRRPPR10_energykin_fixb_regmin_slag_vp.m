% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:27:52
% EndTime: 2019-03-09 16:27:53
% DurationCPUTime: 0.14s
% Computational Cost: add. (614->56), mult. (1405->113), div. (0->0), fcn. (1062->10), ass. (0->48)
t236 = pkin(3) + qJ(5);
t235 = cos(pkin(11));
t213 = sin(pkin(6));
t221 = qJD(1) ^ 2;
t234 = t213 ^ 2 * t221;
t230 = cos(pkin(6)) * qJD(1);
t210 = qJD(2) + t230;
t216 = sin(qJ(3));
t219 = cos(qJ(3));
t217 = sin(qJ(2));
t231 = qJD(1) * t213;
t227 = t217 * t231;
t202 = t216 * t210 + t219 * t227;
t220 = cos(qJ(2));
t226 = t220 * t231;
t205 = -qJD(3) + t226;
t229 = pkin(1) * t230;
t232 = pkin(8) * t226 + t217 * t229;
t197 = t210 * pkin(9) + t232;
t199 = (-pkin(2) * t220 - pkin(9) * t217 - pkin(1)) * t231;
t225 = -t216 * t197 + t219 * t199;
t224 = qJD(4) - t225;
t181 = t202 * pkin(4) + t236 * t205 + t224;
t201 = -t219 * t210 + t216 * t227;
t223 = -pkin(8) * t227 + t220 * t229;
t196 = -t210 * pkin(2) - t223;
t222 = -t202 * qJ(4) + t196;
t183 = t236 * t201 + t222;
t212 = sin(pkin(11));
t177 = t212 * t181 + t235 * t183;
t233 = t219 * t197 + t216 * t199;
t228 = t220 * t234;
t189 = t205 * qJ(4) - t233;
t176 = t235 * t181 - t212 * t183;
t184 = -t201 * pkin(4) + qJD(5) - t189;
t218 = cos(qJ(6));
t215 = sin(qJ(6));
t200 = qJD(6) + t202;
t192 = t212 * t201 - t235 * t205;
t191 = -t235 * t201 - t212 * t205;
t188 = t205 * pkin(3) + t224;
t187 = t201 * pkin(3) + t222;
t186 = -t215 * t191 + t218 * t192;
t185 = t218 * t191 + t215 * t192;
t178 = t191 * pkin(5) + t184;
t175 = -t191 * pkin(10) + t177;
t174 = t202 * pkin(5) - t192 * pkin(10) + t176;
t1 = [t221 / 0.2e1, 0, 0, t217 ^ 2 * t234 / 0.2e1, t217 * t228, t210 * t227, t210 * t226, t210 ^ 2 / 0.2e1, pkin(1) * t228 + t223 * t210, -pkin(1) * t217 * t234 - t232 * t210, t202 ^ 2 / 0.2e1, -t202 * t201, -t202 * t205, t201 * t205, t205 ^ 2 / 0.2e1, t196 * t201 - t225 * t205, t196 * t202 + t233 * t205, t188 * t202 + t189 * t201, -t187 * t201 - t188 * t205, -t187 * t202 + t189 * t205, t187 ^ 2 / 0.2e1 + t189 ^ 2 / 0.2e1 + t188 ^ 2 / 0.2e1, t176 * t202 + t184 * t191, -t177 * t202 + t184 * t192, -t176 * t192 - t177 * t191, t177 ^ 2 / 0.2e1 + t176 ^ 2 / 0.2e1 + t184 ^ 2 / 0.2e1, t186 ^ 2 / 0.2e1, -t186 * t185, t186 * t200, -t185 * t200, t200 ^ 2 / 0.2e1 (t218 * t174 - t215 * t175) * t200 + t178 * t185 -(t215 * t174 + t218 * t175) * t200 + t178 * t186;];
T_reg  = t1;
