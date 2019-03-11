% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:00:00
% EndTime: 2019-03-09 01:00:00
% DurationCPUTime: 0.17s
% Computational Cost: add. (470->51), mult. (1114->112), div. (0->0), fcn. (923->14), ass. (0->51)
t249 = cos(qJ(2));
t262 = qJD(1) * sin(pkin(6));
t227 = qJD(2) * pkin(2) + t249 * t262;
t237 = sin(pkin(7));
t239 = cos(pkin(7));
t261 = qJD(1) * cos(pkin(6));
t251 = t227 * t239 + t237 * t261;
t267 = cos(qJ(4));
t250 = qJD(2) ^ 2;
t265 = t237 ^ 2 * t250;
t235 = t239 * qJD(2) + qJD(3);
t243 = sin(qJ(4));
t244 = sin(qJ(3));
t260 = qJD(2) * t237;
t257 = t244 * t260;
t221 = t243 * t235 + t267 * t257;
t248 = cos(qJ(3));
t255 = t248 * t260;
t232 = -qJD(4) + t255;
t245 = sin(qJ(2));
t226 = pkin(9) * t260 + t245 * t262;
t259 = t248 * t226 + t251 * t244;
t211 = t235 * pkin(10) + t259;
t234 = t239 * t261;
t217 = t234 + (-t227 + (-pkin(3) * t248 - pkin(10) * t244) * qJD(2)) * t237;
t253 = -t243 * t211 + t267 * t217;
t203 = -t232 * pkin(4) - t221 * pkin(11) + t253;
t220 = -t267 * t235 + t243 * t257;
t263 = t267 * t211 + t243 * t217;
t205 = -t220 * pkin(11) + t263;
t242 = sin(qJ(5));
t247 = cos(qJ(5));
t264 = t242 * t203 + t247 * t205;
t256 = (-t237 * t227 + t234) * t260;
t254 = qJD(2) * t262;
t213 = t247 * t220 + t242 * t221;
t252 = t247 * t203 - t242 * t205;
t223 = t244 * t226;
t210 = -t235 * pkin(3) - t251 * t248 + t223;
t206 = t220 * pkin(4) + t210;
t246 = cos(qJ(6));
t241 = sin(qJ(6));
t228 = -qJD(5) + t232;
t214 = -t242 * t220 + t247 * t221;
t212 = qJD(6) + t213;
t208 = t246 * t214 - t241 * t228;
t207 = t241 * t214 + t246 * t228;
t201 = t213 * pkin(5) - t214 * pkin(12) + t206;
t200 = -t228 * pkin(12) + t264;
t199 = t228 * pkin(5) - t252;
t1 = [qJD(1) ^ 2 / 0.2e1, t250 / 0.2e1, t249 * t254, -t245 * t254, t244 ^ 2 * t265 / 0.2e1, t244 * t248 * t265, t235 * t257, t235 * t255, t235 ^ 2 / 0.2e1, -t223 * t235 + (t251 * t235 - t256) * t248, -t259 * t235 + t244 * t256, t221 ^ 2 / 0.2e1, -t221 * t220, -t221 * t232, t220 * t232, t232 ^ 2 / 0.2e1, t210 * t220 - t253 * t232, t210 * t221 + t263 * t232, t214 ^ 2 / 0.2e1, -t214 * t213, -t214 * t228, t213 * t228, t228 ^ 2 / 0.2e1, t206 * t213 - t252 * t228, t206 * t214 + t264 * t228, t208 ^ 2 / 0.2e1, -t208 * t207, t208 * t212, -t207 * t212, t212 ^ 2 / 0.2e1 (-t241 * t200 + t246 * t201) * t212 + t199 * t207 -(t246 * t200 + t241 * t201) * t212 + t199 * t208;];
T_reg  = t1;
