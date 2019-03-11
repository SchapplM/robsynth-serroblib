% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR5
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:46:46
% EndTime: 2019-03-09 13:46:47
% DurationCPUTime: 0.16s
% Computational Cost: add. (661->57), mult. (1816->117), div. (0->0), fcn. (1483->12), ass. (0->52)
t257 = cos(qJ(5));
t234 = sin(pkin(6));
t244 = qJD(1) ^ 2;
t256 = t234 ^ 2 * t244;
t233 = sin(pkin(12));
t235 = cos(pkin(12));
t243 = cos(qJ(2));
t252 = qJD(1) * t234;
t247 = t243 * t252;
t240 = sin(qJ(2));
t248 = t240 * t252;
t223 = -t233 * t248 + t235 * t247;
t222 = qJD(4) - t223;
t251 = cos(pkin(6)) * qJD(1);
t250 = pkin(1) * t251;
t230 = t243 * t250;
t231 = qJD(2) + t251;
t217 = t231 * pkin(2) + t230 + (-pkin(8) - qJ(3)) * t248;
t253 = pkin(8) * t247 + t240 * t250;
t220 = qJ(3) * t247 + t253;
t207 = t233 * t217 + t235 * t220;
t204 = t231 * pkin(9) + t207;
t224 = (t233 * t243 + t235 * t240) * t252;
t225 = qJD(3) + (-pkin(2) * t243 - pkin(1)) * t252;
t211 = -t223 * pkin(3) - t224 * pkin(9) + t225;
t239 = sin(qJ(4));
t242 = cos(qJ(4));
t254 = t242 * t204 + t239 * t211;
t195 = t222 * pkin(10) + t254;
t206 = t235 * t217 - t233 * t220;
t203 = -t231 * pkin(3) - t206;
t215 = t239 * t224 - t242 * t231;
t216 = t242 * t224 + t239 * t231;
t198 = t215 * pkin(4) - t216 * pkin(10) + t203;
t238 = sin(qJ(5));
t255 = t257 * t195 + t238 * t198;
t249 = t243 * t256;
t246 = -t238 * t195 + t257 * t198;
t245 = -t239 * t204 + t242 * t211;
t214 = qJD(5) + t215;
t194 = -t222 * pkin(4) - t245;
t241 = cos(qJ(6));
t237 = sin(qJ(6));
t212 = qJD(6) + t214;
t209 = t257 * t216 + t238 * t222;
t208 = t238 * t216 - t257 * t222;
t200 = -t237 * t208 + t241 * t209;
t199 = t241 * t208 + t237 * t209;
t192 = t208 * pkin(5) + t194;
t191 = -t208 * pkin(11) + t255;
t190 = t214 * pkin(5) - t209 * pkin(11) + t246;
t1 = [t244 / 0.2e1, 0, 0, t240 ^ 2 * t256 / 0.2e1, t240 * t249, t231 * t248, t231 * t247, t231 ^ 2 / 0.2e1, pkin(1) * t249 + (-pkin(8) * t248 + t230) * t231, -pkin(1) * t240 * t256 - t253 * t231, -t206 * t224 + t207 * t223, t207 ^ 2 / 0.2e1 + t206 ^ 2 / 0.2e1 + t225 ^ 2 / 0.2e1, t216 ^ 2 / 0.2e1, -t216 * t215, t216 * t222, -t215 * t222, t222 ^ 2 / 0.2e1, t203 * t215 + t245 * t222, t203 * t216 - t254 * t222, t209 ^ 2 / 0.2e1, -t209 * t208, t209 * t214, -t208 * t214, t214 ^ 2 / 0.2e1, t194 * t208 + t246 * t214, t194 * t209 - t255 * t214, t200 ^ 2 / 0.2e1, -t200 * t199, t200 * t212, -t199 * t212, t212 ^ 2 / 0.2e1 (t241 * t190 - t237 * t191) * t212 + t192 * t199 -(t237 * t190 + t241 * t191) * t212 + t192 * t200;];
T_reg  = t1;
