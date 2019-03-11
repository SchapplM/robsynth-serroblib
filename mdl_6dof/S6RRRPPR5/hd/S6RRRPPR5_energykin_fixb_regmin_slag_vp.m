% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:44:12
% EndTime: 2019-03-09 15:44:12
% DurationCPUTime: 0.20s
% Computational Cost: add. (944->57), mult. (2217->117), div. (0->0), fcn. (1761->12), ass. (0->51)
t255 = cos(qJ(3));
t254 = cos(pkin(12));
t234 = sin(pkin(6));
t242 = qJD(1) ^ 2;
t253 = t234 ^ 2 * t242;
t249 = cos(pkin(6)) * qJD(1);
t230 = qJD(2) + t249;
t238 = sin(qJ(3));
t239 = sin(qJ(2));
t250 = qJD(1) * t234;
t246 = t239 * t250;
t222 = t238 * t230 + t255 * t246;
t241 = cos(qJ(2));
t245 = t241 * t250;
t225 = -qJD(3) + t245;
t248 = pkin(1) * t249;
t251 = pkin(8) * t245 + t239 * t248;
t218 = t230 * pkin(9) + t251;
t220 = (-pkin(2) * t241 - pkin(9) * t239 - pkin(1)) * t250;
t244 = -t238 * t218 + t255 * t220;
t203 = -t225 * pkin(3) - t222 * qJ(4) + t244;
t221 = -t255 * t230 + t238 * t246;
t252 = t255 * t218 + t238 * t220;
t206 = -t221 * qJ(4) + t252;
t233 = sin(pkin(11));
t235 = cos(pkin(11));
t196 = t233 * t203 + t235 * t206;
t194 = -t225 * qJ(5) + t196;
t243 = -pkin(8) * t246 + t241 * t248;
t217 = -t230 * pkin(2) - t243;
t211 = t221 * pkin(3) + qJD(4) + t217;
t212 = t235 * t221 + t233 * t222;
t213 = -t233 * t221 + t235 * t222;
t199 = t212 * pkin(4) - t213 * qJ(5) + t211;
t232 = sin(pkin(12));
t190 = t254 * t194 + t232 * t199;
t247 = t241 * t253;
t189 = -t232 * t194 + t254 * t199;
t195 = t235 * t203 - t233 * t206;
t193 = t225 * pkin(4) + qJD(5) - t195;
t240 = cos(qJ(6));
t237 = sin(qJ(6));
t210 = qJD(6) + t212;
t209 = t254 * t213 - t232 * t225;
t208 = t232 * t213 + t254 * t225;
t201 = -t237 * t208 + t240 * t209;
t200 = t240 * t208 + t237 * t209;
t191 = t208 * pkin(5) + t193;
t188 = -t208 * pkin(10) + t190;
t187 = t212 * pkin(5) - t209 * pkin(10) + t189;
t1 = [t242 / 0.2e1, 0, 0, t239 ^ 2 * t253 / 0.2e1, t239 * t247, t230 * t246, t230 * t245, t230 ^ 2 / 0.2e1, pkin(1) * t247 + t243 * t230, -pkin(1) * t239 * t253 - t251 * t230, t222 ^ 2 / 0.2e1, -t222 * t221, -t222 * t225, t221 * t225, t225 ^ 2 / 0.2e1, t217 * t221 - t244 * t225, t217 * t222 + t252 * t225, -t195 * t213 - t196 * t212, t196 ^ 2 / 0.2e1 + t195 ^ 2 / 0.2e1 + t211 ^ 2 / 0.2e1, t189 * t212 + t193 * t208, -t190 * t212 + t193 * t209, -t189 * t209 - t190 * t208, t190 ^ 2 / 0.2e1 + t189 ^ 2 / 0.2e1 + t193 ^ 2 / 0.2e1, t201 ^ 2 / 0.2e1, -t201 * t200, t201 * t210, -t200 * t210, t210 ^ 2 / 0.2e1 (t240 * t187 - t237 * t188) * t210 + t191 * t200 -(t237 * t187 + t240 * t188) * t210 + t191 * t201;];
T_reg  = t1;
