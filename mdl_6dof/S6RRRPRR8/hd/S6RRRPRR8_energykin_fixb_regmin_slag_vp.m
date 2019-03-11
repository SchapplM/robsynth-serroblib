% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:52:21
% EndTime: 2019-03-09 18:52:21
% DurationCPUTime: 0.19s
% Computational Cost: add. (789->56), mult. (1852->116), div. (0->0), fcn. (1501->12), ass. (0->52)
t256 = cos(qJ(3));
t255 = cos(qJ(5));
t232 = sin(pkin(6));
t241 = qJD(1) ^ 2;
t254 = t232 ^ 2 * t241;
t249 = cos(pkin(6)) * qJD(1);
t229 = qJD(2) + t249;
t237 = sin(qJ(3));
t238 = sin(qJ(2));
t250 = qJD(1) * t232;
t246 = t238 * t250;
t221 = t237 * t229 + t256 * t246;
t240 = cos(qJ(2));
t245 = t240 * t250;
t224 = -qJD(3) + t245;
t248 = pkin(1) * t249;
t251 = pkin(8) * t245 + t238 * t248;
t217 = t229 * pkin(9) + t251;
t219 = (-pkin(2) * t240 - pkin(9) * t238 - pkin(1)) * t250;
t243 = -t237 * t217 + t256 * t219;
t201 = -t224 * pkin(3) - t221 * qJ(4) + t243;
t220 = -t256 * t229 + t237 * t246;
t252 = t256 * t217 + t237 * t219;
t204 = -t220 * qJ(4) + t252;
t231 = sin(pkin(12));
t233 = cos(pkin(12));
t194 = t231 * t201 + t233 * t204;
t192 = -t224 * pkin(10) + t194;
t242 = -pkin(8) * t246 + t240 * t248;
t216 = -t229 * pkin(2) - t242;
t210 = t220 * pkin(3) + qJD(4) + t216;
t211 = -t233 * t220 - t231 * t221;
t212 = -t231 * t220 + t233 * t221;
t197 = -t211 * pkin(4) - t212 * pkin(10) + t210;
t236 = sin(qJ(5));
t253 = t255 * t192 + t236 * t197;
t247 = t240 * t254;
t244 = -t236 * t192 + t255 * t197;
t193 = t233 * t201 - t231 * t204;
t209 = qJD(5) - t211;
t191 = t224 * pkin(4) - t193;
t239 = cos(qJ(6));
t235 = sin(qJ(6));
t208 = qJD(6) + t209;
t207 = t255 * t212 - t236 * t224;
t206 = t236 * t212 + t255 * t224;
t199 = -t235 * t206 + t239 * t207;
t198 = t239 * t206 + t235 * t207;
t189 = t206 * pkin(5) + t191;
t188 = -t206 * pkin(11) + t253;
t187 = t209 * pkin(5) - t207 * pkin(11) + t244;
t1 = [t241 / 0.2e1, 0, 0, t238 ^ 2 * t254 / 0.2e1, t238 * t247, t229 * t246, t229 * t245, t229 ^ 2 / 0.2e1, pkin(1) * t247 + t242 * t229, -pkin(1) * t238 * t254 - t251 * t229, t221 ^ 2 / 0.2e1, -t221 * t220, -t221 * t224, t220 * t224, t224 ^ 2 / 0.2e1, t216 * t220 - t243 * t224, t216 * t221 + t252 * t224, -t193 * t212 + t194 * t211, t194 ^ 2 / 0.2e1 + t193 ^ 2 / 0.2e1 + t210 ^ 2 / 0.2e1, t207 ^ 2 / 0.2e1, -t207 * t206, t207 * t209, -t206 * t209, t209 ^ 2 / 0.2e1, t191 * t206 + t244 * t209, t191 * t207 - t253 * t209, t199 ^ 2 / 0.2e1, -t199 * t198, t199 * t208, -t198 * t208, t208 ^ 2 / 0.2e1 (t239 * t187 - t235 * t188) * t208 + t189 * t198 -(t235 * t187 + t239 * t188) * t208 + t189 * t199;];
T_reg  = t1;
