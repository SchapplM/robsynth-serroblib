% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:11:06
% EndTime: 2019-03-09 17:11:06
% DurationCPUTime: 0.18s
% Computational Cost: add. (830->53), mult. (1943->109), div. (0->0), fcn. (1527->10), ass. (0->47)
t245 = cos(qJ(3));
t223 = sin(pkin(6));
t231 = qJD(1) ^ 2;
t244 = t223 ^ 2 * t231;
t239 = cos(pkin(6)) * qJD(1);
t220 = qJD(2) + t239;
t227 = sin(qJ(3));
t228 = sin(qJ(2));
t240 = qJD(1) * t223;
t236 = t228 * t240;
t212 = t227 * t220 + t245 * t236;
t230 = cos(qJ(2));
t235 = t230 * t240;
t215 = -qJD(3) + t235;
t238 = pkin(1) * t239;
t241 = pkin(8) * t235 + t228 * t238;
t208 = pkin(9) * t220 + t241;
t210 = (-pkin(2) * t230 - pkin(9) * t228 - pkin(1)) * t240;
t234 = -t208 * t227 + t245 * t210;
t194 = -pkin(3) * t215 - qJ(4) * t212 + t234;
t211 = -t245 * t220 + t227 * t236;
t242 = t245 * t208 + t227 * t210;
t197 = -qJ(4) * t211 + t242;
t222 = sin(pkin(11));
t224 = cos(pkin(11));
t190 = t222 * t194 + t224 * t197;
t188 = -pkin(10) * t215 + t190;
t232 = -pkin(8) * t236 + t230 * t238;
t207 = -t220 * pkin(2) - t232;
t201 = t211 * pkin(3) + qJD(4) + t207;
t202 = -t224 * t211 - t212 * t222;
t203 = -t211 * t222 + t212 * t224;
t192 = -t202 * pkin(4) - t203 * pkin(10) + t201;
t226 = sin(qJ(5));
t229 = cos(qJ(5));
t243 = t229 * t188 + t226 * t192;
t237 = t230 * t244;
t189 = t194 * t224 - t222 * t197;
t233 = -t188 * t226 + t192 * t229;
t187 = pkin(4) * t215 - t189;
t200 = qJD(5) - t202;
t199 = t203 * t229 - t215 * t226;
t198 = t203 * t226 + t229 * t215;
t185 = pkin(5) * t198 - qJ(6) * t199 + t187;
t184 = qJ(6) * t200 + t243;
t183 = -pkin(5) * t200 + qJD(6) - t233;
t1 = [t231 / 0.2e1, 0, 0, t228 ^ 2 * t244 / 0.2e1, t228 * t237, t220 * t236, t220 * t235, t220 ^ 2 / 0.2e1, pkin(1) * t237 + t232 * t220, -pkin(1) * t228 * t244 - t241 * t220, t212 ^ 2 / 0.2e1, -t212 * t211, -t212 * t215, t211 * t215, t215 ^ 2 / 0.2e1, t207 * t211 - t234 * t215, t207 * t212 + t242 * t215, -t189 * t203 + t190 * t202, t190 ^ 2 / 0.2e1 + t189 ^ 2 / 0.2e1 + t201 ^ 2 / 0.2e1, t199 ^ 2 / 0.2e1, -t199 * t198, t199 * t200, -t198 * t200, t200 ^ 2 / 0.2e1, t187 * t198 + t233 * t200, t187 * t199 - t243 * t200, -t183 * t200 + t185 * t198, t183 * t199 - t184 * t198, t184 * t200 - t185 * t199, t184 ^ 2 / 0.2e1 + t185 ^ 2 / 0.2e1 + t183 ^ 2 / 0.2e1;];
T_reg  = t1;
