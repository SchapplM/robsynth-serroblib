% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:08:30
% EndTime: 2019-03-09 23:08:30
% DurationCPUTime: 0.14s
% Computational Cost: add. (645->55), mult. (1492->112), div. (0->0), fcn. (1169->10), ass. (0->49)
t236 = pkin(4) + pkin(11);
t235 = cos(qJ(3));
t210 = sin(pkin(6));
t219 = qJD(1) ^ 2;
t234 = t210 ^ 2 * t219;
t229 = cos(pkin(6)) * qJD(1);
t208 = qJD(2) + t229;
t214 = sin(qJ(3));
t215 = sin(qJ(2));
t230 = qJD(1) * t210;
t226 = t215 * t230;
t199 = t214 * t208 + t235 * t226;
t218 = cos(qJ(2));
t225 = t218 * t230;
t203 = -qJD(3) + t225;
t228 = pkin(1) * t229;
t231 = pkin(8) * t225 + t215 * t228;
t195 = t208 * pkin(9) + t231;
t197 = (-pkin(2) * t218 - pkin(9) * t215 - pkin(1)) * t230;
t223 = -t214 * t195 + t235 * t197;
t180 = -t203 * pkin(3) - t199 * pkin(10) + t223;
t198 = -t235 * t208 + t214 * t226;
t232 = t235 * t195 + t214 * t197;
t183 = -t198 * pkin(10) + t232;
t213 = sin(qJ(4));
t217 = cos(qJ(4));
t233 = t213 * t180 + t217 * t183;
t227 = t218 * t234;
t224 = t217 * t180 - t213 * t183;
t200 = -qJD(4) + t203;
t177 = t200 * qJ(5) - t233;
t222 = qJD(5) - t224;
t189 = -t213 * t198 + t217 * t199;
t221 = -pkin(8) * t226 + t218 * t228;
t194 = -t208 * pkin(2) - t221;
t190 = t198 * pkin(3) + t194;
t220 = -t189 * qJ(5) + t190;
t216 = cos(qJ(6));
t212 = sin(qJ(6));
t188 = t217 * t198 + t213 * t199;
t187 = qJD(6) + t189;
t185 = t212 * t188 - t216 * t200;
t184 = -t216 * t188 - t212 * t200;
t178 = t188 * pkin(4) + t220;
t176 = t200 * pkin(4) + t222;
t175 = t236 * t188 + t220;
t174 = -t188 * pkin(5) - t177;
t173 = t189 * pkin(5) + t236 * t200 + t222;
t1 = [t219 / 0.2e1, 0, 0, t215 ^ 2 * t234 / 0.2e1, t215 * t227, t208 * t226, t208 * t225, t208 ^ 2 / 0.2e1, pkin(1) * t227 + t221 * t208, -pkin(1) * t215 * t234 - t231 * t208, t199 ^ 2 / 0.2e1, -t199 * t198, -t199 * t203, t198 * t203, t203 ^ 2 / 0.2e1, t194 * t198 - t223 * t203, t194 * t199 + t232 * t203, t189 ^ 2 / 0.2e1, -t189 * t188, -t189 * t200, t188 * t200, t200 ^ 2 / 0.2e1, t190 * t188 - t224 * t200, t190 * t189 + t233 * t200, t176 * t189 + t177 * t188, -t176 * t200 - t178 * t188, t177 * t200 - t178 * t189, t178 ^ 2 / 0.2e1 + t177 ^ 2 / 0.2e1 + t176 ^ 2 / 0.2e1, t185 ^ 2 / 0.2e1, -t185 * t184, t185 * t187, -t184 * t187, t187 ^ 2 / 0.2e1 (t216 * t173 - t212 * t175) * t187 + t174 * t184 -(t212 * t173 + t216 * t175) * t187 + t174 * t185;];
T_reg  = t1;
