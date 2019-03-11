% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:28:17
% EndTime: 2019-03-08 23:28:18
% DurationCPUTime: 0.20s
% Computational Cost: add. (481->50), mult. (1170->109), div. (0->0), fcn. (956->14), ass. (0->50)
t252 = cos(qJ(2));
t264 = qJD(1) * sin(pkin(6));
t232 = qJD(2) * pkin(2) + t252 * t264;
t241 = sin(pkin(7));
t244 = cos(pkin(7));
t263 = qJD(1) * cos(pkin(6));
t254 = t232 * t244 + t241 * t263;
t268 = cos(qJ(4));
t253 = qJD(2) ^ 2;
t266 = t241 ^ 2 * t253;
t238 = t244 * qJD(2) + qJD(3);
t247 = sin(qJ(4));
t248 = sin(qJ(3));
t262 = qJD(2) * t241;
t259 = t248 * t262;
t226 = t247 * t238 + t268 * t259;
t251 = cos(qJ(3));
t257 = t251 * t262;
t235 = -qJD(4) + t257;
t249 = sin(qJ(2));
t230 = pkin(9) * t262 + t249 * t264;
t261 = t251 * t230 + t254 * t248;
t216 = t238 * pkin(10) + t261;
t237 = t244 * t263;
t222 = t237 + (-t232 + (-pkin(3) * t251 - pkin(10) * t248) * qJD(2)) * t241;
t255 = -t247 * t216 + t268 * t222;
t208 = -t235 * pkin(4) - t226 * qJ(5) + t255;
t225 = -t268 * t238 + t247 * t259;
t265 = t268 * t216 + t247 * t222;
t210 = -t225 * qJ(5) + t265;
t240 = sin(pkin(13));
t243 = cos(pkin(13));
t205 = t240 * t208 + t243 * t210;
t258 = (-t241 * t232 + t237) * t262;
t256 = qJD(2) * t264;
t218 = -t243 * t225 - t240 * t226;
t204 = t243 * t208 - t240 * t210;
t228 = t248 * t230;
t215 = -t238 * pkin(3) - t254 * t251 + t228;
t211 = t225 * pkin(4) + qJD(5) + t215;
t250 = cos(qJ(6));
t246 = sin(qJ(6));
t219 = -t240 * t225 + t243 * t226;
t217 = qJD(6) - t218;
t213 = t250 * t219 - t246 * t235;
t212 = t246 * t219 + t250 * t235;
t206 = -t218 * pkin(5) - t219 * pkin(11) + t211;
t203 = -t235 * pkin(11) + t205;
t202 = t235 * pkin(5) - t204;
t1 = [qJD(1) ^ 2 / 0.2e1, t253 / 0.2e1, t252 * t256, -t249 * t256, t248 ^ 2 * t266 / 0.2e1, t248 * t251 * t266, t238 * t259, t238 * t257, t238 ^ 2 / 0.2e1, -t228 * t238 + (t254 * t238 - t258) * t251, -t261 * t238 + t248 * t258, t226 ^ 2 / 0.2e1, -t226 * t225, -t226 * t235, t225 * t235, t235 ^ 2 / 0.2e1, t215 * t225 - t255 * t235, t215 * t226 + t265 * t235, -t204 * t219 + t205 * t218, t205 ^ 2 / 0.2e1 + t204 ^ 2 / 0.2e1 + t211 ^ 2 / 0.2e1, t213 ^ 2 / 0.2e1, -t213 * t212, t213 * t217, -t212 * t217, t217 ^ 2 / 0.2e1 (-t246 * t203 + t250 * t206) * t217 + t202 * t212 -(t250 * t203 + t246 * t206) * t217 + t202 * t213;];
T_reg  = t1;
