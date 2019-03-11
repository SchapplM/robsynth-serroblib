% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR7
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
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:39:25
% EndTime: 2019-03-09 18:39:25
% DurationCPUTime: 0.18s
% Computational Cost: add. (857->56), mult. (2024->116), div. (0->0), fcn. (1651->12), ass. (0->52)
t253 = cos(qJ(3));
t229 = sin(pkin(6));
t239 = qJD(1) ^ 2;
t252 = t229 ^ 2 * t239;
t247 = cos(pkin(6)) * qJD(1);
t226 = qJD(2) + t247;
t234 = sin(qJ(3));
t235 = sin(qJ(2));
t248 = qJD(1) * t229;
t244 = t235 * t248;
t216 = t234 * t226 + t253 * t244;
t238 = cos(qJ(2));
t243 = t238 * t248;
t221 = -qJD(3) + t243;
t246 = pkin(1) * t247;
t249 = pkin(8) * t243 + t235 * t246;
t213 = t226 * pkin(9) + t249;
t214 = (-pkin(2) * t238 - pkin(9) * t235 - pkin(1)) * t248;
t242 = -t234 * t213 + t253 * t214;
t202 = -t221 * pkin(3) - t216 * qJ(4) + t242;
t215 = -t253 * t226 + t234 * t244;
t250 = t253 * t213 + t234 * t214;
t204 = -t215 * qJ(4) + t250;
t228 = sin(pkin(12));
t230 = cos(pkin(12));
t192 = t230 * t202 - t228 * t204;
t208 = -t228 * t215 + t230 * t216;
t189 = -t221 * pkin(4) - t208 * pkin(10) + t192;
t193 = t228 * t202 + t230 * t204;
t207 = -t230 * t215 - t228 * t216;
t191 = t207 * pkin(10) + t193;
t233 = sin(qJ(5));
t237 = cos(qJ(5));
t251 = t233 * t189 + t237 * t191;
t245 = t238 * t252;
t197 = -t237 * t207 + t233 * t208;
t241 = t237 * t189 - t233 * t191;
t240 = -pkin(8) * t244 + t238 * t246;
t212 = -t226 * pkin(2) - t240;
t206 = t215 * pkin(3) + qJD(4) + t212;
t199 = -t207 * pkin(4) + t206;
t236 = cos(qJ(6));
t232 = sin(qJ(6));
t218 = -qJD(5) + t221;
t198 = t233 * t207 + t237 * t208;
t196 = qJD(6) + t197;
t195 = t236 * t198 - t232 * t218;
t194 = t232 * t198 + t236 * t218;
t187 = t197 * pkin(5) - t198 * pkin(11) + t199;
t186 = -t218 * pkin(11) + t251;
t185 = t218 * pkin(5) - t241;
t1 = [t239 / 0.2e1, 0, 0, t235 ^ 2 * t252 / 0.2e1, t235 * t245, t226 * t244, t226 * t243, t226 ^ 2 / 0.2e1, pkin(1) * t245 + t240 * t226, -pkin(1) * t235 * t252 - t249 * t226, t216 ^ 2 / 0.2e1, -t216 * t215, -t216 * t221, t215 * t221, t221 ^ 2 / 0.2e1, t212 * t215 - t242 * t221, t212 * t216 + t250 * t221, -t192 * t208 + t193 * t207, t193 ^ 2 / 0.2e1 + t192 ^ 2 / 0.2e1 + t206 ^ 2 / 0.2e1, t198 ^ 2 / 0.2e1, -t198 * t197, -t198 * t218, t197 * t218, t218 ^ 2 / 0.2e1, t199 * t197 - t241 * t218, t199 * t198 + t218 * t251, t195 ^ 2 / 0.2e1, -t195 * t194, t195 * t196, -t194 * t196, t196 ^ 2 / 0.2e1 (-t232 * t186 + t236 * t187) * t196 + t185 * t194 -(t236 * t186 + t232 * t187) * t196 + t185 * t195;];
T_reg  = t1;
