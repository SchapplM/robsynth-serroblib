% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRPRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:37
% EndTime: 2019-03-08 22:40:37
% DurationCPUTime: 0.17s
% Computational Cost: add. (290->53), mult. (717->112), div. (0->0), fcn. (552->12), ass. (0->49)
t216 = sin(pkin(7));
t239 = qJD(1) * cos(pkin(6));
t227 = cos(qJ(2));
t240 = qJD(1) * sin(pkin(6));
t208 = qJD(2) * pkin(2) + t227 * t240;
t218 = cos(pkin(7));
t243 = t208 * t218;
t229 = t216 * t239 + t243;
t245 = -pkin(3) - pkin(10);
t222 = sin(qJ(3));
t244 = qJ(4) * t222;
t228 = qJD(2) ^ 2;
t242 = t216 ^ 2 * t228;
t214 = t218 * qJD(2) + qJD(3);
t226 = cos(qJ(3));
t223 = sin(qJ(2));
t238 = qJD(2) * t216;
t206 = pkin(9) * t238 + t223 * t240;
t204 = t222 * t206;
t236 = qJD(4) + t204;
t237 = qJD(2) * t222;
t190 = -t226 * t243 + (pkin(4) * t237 - t226 * t239) * t216 + t245 * t214 + t236;
t213 = t218 * t239;
t195 = t213 + (-t208 + (t245 * t226 - t244) * qJD(2)) * t216;
t221 = sin(qJ(5));
t225 = cos(qJ(5));
t241 = t221 * t190 + t225 * t195;
t235 = t226 * t206 + t229 * t222;
t233 = t226 * t238;
t232 = t216 * t237;
t231 = qJD(2) * t240;
t193 = -t214 * qJ(4) - t235;
t191 = pkin(4) * t233 - t193;
t230 = t225 * t190 - t221 * t195;
t201 = t221 * t214 + t225 * t233;
t224 = cos(qJ(6));
t220 = sin(qJ(6));
t210 = qJD(5) + t232;
t202 = t225 * t214 - t221 * t233;
t200 = qJD(6) + t201;
t199 = -t216 * t208 + t213;
t198 = t224 * t202 + t220 * t210;
t197 = t220 * t202 - t224 * t210;
t196 = t213 + (-t208 + (-pkin(3) * t226 - t244) * qJD(2)) * t216;
t192 = -t214 * pkin(3) - t229 * t226 + t236;
t188 = t201 * pkin(5) - t202 * pkin(11) + t191;
t187 = t210 * pkin(11) + t241;
t186 = -t210 * pkin(5) - t230;
t1 = [qJD(1) ^ 2 / 0.2e1, t228 / 0.2e1, t227 * t231, -t223 * t231, t222 ^ 2 * t242 / 0.2e1, t222 * t226 * t242, t214 * t232, t214 * t233, t214 ^ 2 / 0.2e1, -t204 * t214 + (-t199 * t238 + t229 * t214) * t226, t199 * t232 - t235 * t214 (t192 * t222 - t193 * t226) * t238, t192 * t214 + t196 * t233, -t193 * t214 - t196 * t232, t196 ^ 2 / 0.2e1 + t193 ^ 2 / 0.2e1 + t192 ^ 2 / 0.2e1, t202 ^ 2 / 0.2e1, -t202 * t201, t202 * t210, -t201 * t210, t210 ^ 2 / 0.2e1, t191 * t201 + t230 * t210, t191 * t202 - t241 * t210, t198 ^ 2 / 0.2e1, -t198 * t197, t198 * t200, -t197 * t200, t200 ^ 2 / 0.2e1 (-t220 * t187 + t224 * t188) * t200 + t186 * t197 -(t224 * t187 + t220 * t188) * t200 + t186 * t198;];
T_reg  = t1;
