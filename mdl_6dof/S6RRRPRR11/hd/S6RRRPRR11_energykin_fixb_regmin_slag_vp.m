% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:34:14
% EndTime: 2019-03-09 19:34:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (518->55), mult. (1194->112), div. (0->0), fcn. (916->10), ass. (0->47)
t218 = sin(pkin(6));
t228 = qJD(1) ^ 2;
t242 = t218 ^ 2 * t228;
t236 = cos(pkin(6)) * qJD(1);
t216 = qJD(2) + t236;
t222 = sin(qJ(3));
t226 = cos(qJ(3));
t223 = sin(qJ(2));
t237 = qJD(1) * t218;
t233 = t223 * t237;
t204 = t222 * t216 + t226 * t233;
t227 = cos(qJ(2));
t232 = t227 * t237;
t209 = -qJD(3) + t232;
t235 = pkin(1) * t236;
t238 = pkin(8) * t232 + t223 * t235;
t198 = t216 * pkin(9) + t238;
t200 = (-pkin(2) * t227 - pkin(9) * t223 - pkin(1)) * t237;
t231 = -t222 * t198 + t226 * t200;
t230 = qJD(4) - t231;
t182 = -t204 * pkin(10) + (pkin(3) + pkin(4)) * t209 + t230;
t240 = t226 * t198 + t222 * t200;
t188 = -t209 * qJ(4) + t240;
t203 = -t226 * t216 + t222 * t233;
t185 = t203 * pkin(10) + t188;
t221 = sin(qJ(5));
t225 = cos(qJ(5));
t241 = t221 * t182 + t225 * t185;
t239 = -pkin(8) * t233 + t227 * t235;
t234 = t227 * t242;
t197 = -t216 * pkin(2) - t239;
t192 = -t225 * t203 + t221 * t204;
t186 = t203 * pkin(3) - t204 * qJ(4) + t197;
t229 = t225 * t182 - t221 * t185;
t183 = -t203 * pkin(4) - t186;
t224 = cos(qJ(6));
t220 = sin(qJ(6));
t207 = qJD(5) + t209;
t193 = t221 * t203 + t225 * t204;
t191 = qJD(6) + t192;
t190 = t224 * t193 + t220 * t207;
t189 = t220 * t193 - t224 * t207;
t187 = t209 * pkin(3) + t230;
t180 = t192 * pkin(5) - t193 * pkin(11) + t183;
t179 = t207 * pkin(11) + t241;
t178 = -t207 * pkin(5) - t229;
t1 = [t228 / 0.2e1, 0, 0, t223 ^ 2 * t242 / 0.2e1, t223 * t234, t216 * t233, t216 * t232, t216 ^ 2 / 0.2e1, pkin(1) * t234 + t239 * t216, -pkin(1) * t223 * t242 - t238 * t216, t204 ^ 2 / 0.2e1, -t204 * t203, -t204 * t209, t203 * t209, t209 ^ 2 / 0.2e1, t197 * t203 - t231 * t209, t197 * t204 + t240 * t209, t186 * t203 + t187 * t209, t187 * t204 - t188 * t203, -t186 * t204 - t188 * t209, t188 ^ 2 / 0.2e1 + t186 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1, t193 ^ 2 / 0.2e1, -t193 * t192, t193 * t207, -t192 * t207, t207 ^ 2 / 0.2e1, t183 * t192 + t229 * t207, t183 * t193 - t241 * t207, t190 ^ 2 / 0.2e1, -t190 * t189, t190 * t191, -t189 * t191, t191 ^ 2 / 0.2e1 (-t220 * t179 + t224 * t180) * t191 + t178 * t189 -(t224 * t179 + t220 * t180) * t191 + t178 * t190;];
T_reg  = t1;
