% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:35:10
% EndTime: 2019-03-09 10:35:10
% DurationCPUTime: 0.17s
% Computational Cost: add. (797->58), mult. (2183->118), div. (0->0), fcn. (1744->12), ass. (0->51)
t255 = cos(pkin(12));
t235 = sin(pkin(6));
t244 = qJD(1) ^ 2;
t254 = t235 ^ 2 * t244;
t234 = sin(pkin(11));
t236 = cos(pkin(11));
t243 = cos(qJ(2));
t251 = qJD(1) * t235;
t246 = t243 * t251;
t240 = sin(qJ(2));
t247 = t240 * t251;
t223 = -t234 * t247 + t236 * t246;
t222 = qJD(4) - t223;
t250 = cos(pkin(6)) * qJD(1);
t249 = pkin(1) * t250;
t230 = t243 * t249;
t231 = qJD(2) + t250;
t217 = t231 * pkin(2) + t230 + (-pkin(8) - qJ(3)) * t247;
t252 = pkin(8) * t246 + t240 * t249;
t221 = qJ(3) * t246 + t252;
t208 = t234 * t217 + t236 * t221;
t205 = pkin(9) * t231 + t208;
t224 = (t234 * t243 + t236 * t240) * t251;
t225 = qJD(3) + (-pkin(2) * t243 - pkin(1)) * t251;
t212 = -t223 * pkin(3) - t224 * pkin(9) + t225;
t239 = sin(qJ(4));
t242 = cos(qJ(4));
t253 = t242 * t205 + t239 * t212;
t196 = qJ(5) * t222 + t253;
t207 = t217 * t236 - t234 * t221;
t204 = -pkin(3) * t231 - t207;
t215 = t224 * t239 - t242 * t231;
t216 = t224 * t242 + t231 * t239;
t199 = pkin(4) * t215 - qJ(5) * t216 + t204;
t233 = sin(pkin(12));
t192 = t255 * t196 + t233 * t199;
t248 = t243 * t254;
t191 = -t196 * t233 + t255 * t199;
t245 = -t239 * t205 + t212 * t242;
t195 = -pkin(4) * t222 + qJD(5) - t245;
t241 = cos(qJ(6));
t238 = sin(qJ(6));
t214 = qJD(6) + t215;
t210 = t255 * t216 + t233 * t222;
t209 = t233 * t216 - t255 * t222;
t201 = -t209 * t238 + t210 * t241;
t200 = t241 * t209 + t210 * t238;
t193 = pkin(5) * t209 + t195;
t190 = -pkin(10) * t209 + t192;
t189 = pkin(5) * t215 - pkin(10) * t210 + t191;
t1 = [t244 / 0.2e1, 0, 0, t240 ^ 2 * t254 / 0.2e1, t240 * t248, t231 * t247, t231 * t246, t231 ^ 2 / 0.2e1, pkin(1) * t248 + (-pkin(8) * t247 + t230) * t231, -pkin(1) * t240 * t254 - t252 * t231, -t207 * t224 + t208 * t223, t208 ^ 2 / 0.2e1 + t207 ^ 2 / 0.2e1 + t225 ^ 2 / 0.2e1, t216 ^ 2 / 0.2e1, -t216 * t215, t216 * t222, -t215 * t222, t222 ^ 2 / 0.2e1, t204 * t215 + t245 * t222, t204 * t216 - t253 * t222, t191 * t215 + t195 * t209, -t192 * t215 + t195 * t210, -t191 * t210 - t192 * t209, t192 ^ 2 / 0.2e1 + t191 ^ 2 / 0.2e1 + t195 ^ 2 / 0.2e1, t201 ^ 2 / 0.2e1, -t201 * t200, t201 * t214, -t200 * t214, t214 ^ 2 / 0.2e1 (t189 * t241 - t190 * t238) * t214 + t193 * t200 -(t189 * t238 + t190 * t241) * t214 + t193 * t201;];
T_reg  = t1;
