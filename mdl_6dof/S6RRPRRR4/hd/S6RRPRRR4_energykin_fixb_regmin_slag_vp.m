% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:34:48
% EndTime: 2019-03-09 13:34:48
% DurationCPUTime: 0.17s
% Computational Cost: add. (678->57), mult. (1859->117), div. (0->0), fcn. (1518->12), ass. (0->52)
t259 = cos(qJ(4));
t236 = sin(pkin(6));
t246 = qJD(1) ^ 2;
t258 = t236 ^ 2 * t246;
t235 = sin(pkin(12));
t237 = cos(pkin(12));
t242 = sin(qJ(2));
t245 = cos(qJ(2));
t254 = qJD(1) * t236;
t226 = (t235 * t245 + t237 * t242) * t254;
t253 = cos(pkin(6)) * qJD(1);
t233 = qJD(2) + t253;
t241 = sin(qJ(4));
t217 = t259 * t226 + t241 * t233;
t249 = t245 * t254;
t250 = t242 * t254;
t225 = -t235 * t250 + t237 * t249;
t224 = qJD(4) - t225;
t252 = pkin(1) * t253;
t232 = t245 * t252;
t218 = t233 * pkin(2) + t232 + (-pkin(8) - qJ(3)) * t250;
t255 = pkin(8) * t249 + t242 * t252;
t222 = qJ(3) * t249 + t255;
t210 = t235 * t218 + t237 * t222;
t208 = t233 * pkin(9) + t210;
t227 = qJD(3) + (-pkin(2) * t245 - pkin(1)) * t254;
t213 = -t225 * pkin(3) - t226 * pkin(9) + t227;
t248 = -t241 * t208 + t259 * t213;
t197 = t224 * pkin(4) - t217 * pkin(10) + t248;
t216 = t241 * t226 - t259 * t233;
t256 = t259 * t208 + t241 * t213;
t199 = -t216 * pkin(10) + t256;
t240 = sin(qJ(5));
t244 = cos(qJ(5));
t257 = t240 * t197 + t244 * t199;
t251 = t245 * t258;
t205 = t244 * t216 + t240 * t217;
t209 = t237 * t218 - t235 * t222;
t247 = t244 * t197 - t240 * t199;
t207 = -t233 * pkin(3) - t209;
t200 = t216 * pkin(4) + t207;
t243 = cos(qJ(6));
t239 = sin(qJ(6));
t223 = qJD(5) + t224;
t206 = -t240 * t216 + t244 * t217;
t204 = qJD(6) + t205;
t202 = t243 * t206 + t239 * t223;
t201 = t239 * t206 - t243 * t223;
t195 = t205 * pkin(5) - t206 * pkin(11) + t200;
t194 = t223 * pkin(11) + t257;
t193 = -t223 * pkin(5) - t247;
t1 = [t246 / 0.2e1, 0, 0, t242 ^ 2 * t258 / 0.2e1, t242 * t251, t233 * t250, t233 * t249, t233 ^ 2 / 0.2e1, pkin(1) * t251 + (-pkin(8) * t250 + t232) * t233, -pkin(1) * t242 * t258 - t255 * t233, -t209 * t226 + t210 * t225, t210 ^ 2 / 0.2e1 + t209 ^ 2 / 0.2e1 + t227 ^ 2 / 0.2e1, t217 ^ 2 / 0.2e1, -t217 * t216, t217 * t224, -t216 * t224, t224 ^ 2 / 0.2e1, t207 * t216 + t248 * t224, t207 * t217 - t256 * t224, t206 ^ 2 / 0.2e1, -t206 * t205, t206 * t223, -t205 * t223, t223 ^ 2 / 0.2e1, t200 * t205 + t247 * t223, t200 * t206 - t257 * t223, t202 ^ 2 / 0.2e1, -t202 * t201, t202 * t204, -t201 * t204, t204 ^ 2 / 0.2e1 (-t239 * t194 + t243 * t195) * t204 + t193 * t201 -(t243 * t194 + t239 * t195) * t204 + t193 * t202;];
T_reg  = t1;
