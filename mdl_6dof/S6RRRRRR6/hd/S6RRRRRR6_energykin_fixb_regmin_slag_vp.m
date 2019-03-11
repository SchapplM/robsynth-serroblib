% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:21:34
% EndTime: 2019-03-10 04:21:34
% DurationCPUTime: 0.17s
% Computational Cost: add. (817->57), mult. (1833->119), div. (0->0), fcn. (1505->12), ass. (0->53)
t255 = cos(qJ(3));
t254 = cos(qJ(5));
t228 = sin(pkin(6));
t238 = qJD(1) ^ 2;
t253 = t228 ^ 2 * t238;
t237 = cos(qJ(2));
t248 = qJD(1) * t228;
t243 = t237 * t248;
t221 = -qJD(3) + t243;
t218 = -qJD(4) + t221;
t247 = cos(pkin(6)) * qJD(1);
t226 = qJD(2) + t247;
t233 = sin(qJ(3));
t234 = sin(qJ(2));
t244 = t234 * t248;
t216 = t233 * t226 + t255 * t244;
t246 = pkin(1) * t247;
t249 = pkin(8) * t243 + t234 * t246;
t212 = t226 * pkin(9) + t249;
t214 = (-pkin(2) * t237 - pkin(9) * t234 - pkin(1)) * t248;
t240 = -t233 * t212 + t255 * t214;
t196 = -t221 * pkin(3) - t216 * pkin(10) + t240;
t215 = -t255 * t226 + t233 * t244;
t250 = t255 * t212 + t233 * t214;
t199 = -t215 * pkin(10) + t250;
t232 = sin(qJ(4));
t236 = cos(qJ(4));
t251 = t232 * t196 + t236 * t199;
t189 = -t218 * pkin(11) + t251;
t205 = t236 * t215 + t232 * t216;
t206 = -t232 * t215 + t236 * t216;
t239 = -pkin(8) * t244 + t237 * t246;
t211 = -t226 * pkin(2) - t239;
t207 = t215 * pkin(3) + t211;
t192 = t205 * pkin(4) - t206 * pkin(11) + t207;
t231 = sin(qJ(5));
t252 = t254 * t189 + t231 * t192;
t245 = t237 * t253;
t242 = -t231 * t189 + t254 * t192;
t241 = t236 * t196 - t232 * t199;
t204 = qJD(5) + t205;
t188 = t218 * pkin(4) - t241;
t235 = cos(qJ(6));
t230 = sin(qJ(6));
t203 = qJD(6) + t204;
t202 = t254 * t206 - t231 * t218;
t201 = t231 * t206 + t254 * t218;
t194 = -t230 * t201 + t235 * t202;
t193 = t235 * t201 + t230 * t202;
t186 = t201 * pkin(5) + t188;
t185 = -t201 * pkin(12) + t252;
t184 = t204 * pkin(5) - t202 * pkin(12) + t242;
t1 = [t238 / 0.2e1, 0, 0, t234 ^ 2 * t253 / 0.2e1, t234 * t245, t226 * t244, t226 * t243, t226 ^ 2 / 0.2e1, pkin(1) * t245 + t239 * t226, -pkin(1) * t234 * t253 - t249 * t226, t216 ^ 2 / 0.2e1, -t216 * t215, -t216 * t221, t215 * t221, t221 ^ 2 / 0.2e1, t211 * t215 - t240 * t221, t211 * t216 + t250 * t221, t206 ^ 2 / 0.2e1, -t206 * t205, -t206 * t218, t205 * t218, t218 ^ 2 / 0.2e1, t207 * t205 - t241 * t218, t207 * t206 + t251 * t218, t202 ^ 2 / 0.2e1, -t202 * t201, t202 * t204, -t201 * t204, t204 ^ 2 / 0.2e1, t188 * t201 + t242 * t204, t188 * t202 - t252 * t204, t194 ^ 2 / 0.2e1, -t194 * t193, t194 * t203, -t193 * t203, t203 ^ 2 / 0.2e1 (t235 * t184 - t230 * t185) * t203 + t186 * t193 -(t230 * t184 + t235 * t185) * t203 + t186 * t194;];
T_reg  = t1;
