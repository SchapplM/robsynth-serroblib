% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:55:43
% EndTime: 2019-03-09 22:55:43
% DurationCPUTime: 0.18s
% Computational Cost: add. (970->58), mult. (2198->120), div. (0->0), fcn. (1765->12), ass. (0->52)
t253 = cos(qJ(3));
t252 = cos(pkin(12));
t229 = sin(pkin(6));
t238 = qJD(1) ^ 2;
t251 = t229 ^ 2 * t238;
t237 = cos(qJ(2));
t247 = qJD(1) * t229;
t242 = t237 * t247;
t221 = -qJD(3) + t242;
t218 = -qJD(4) + t221;
t246 = cos(pkin(6)) * qJD(1);
t226 = qJD(2) + t246;
t233 = sin(qJ(3));
t234 = sin(qJ(2));
t243 = t234 * t247;
t216 = t233 * t226 + t253 * t243;
t245 = pkin(1) * t246;
t248 = pkin(8) * t242 + t234 * t245;
t212 = t226 * pkin(9) + t248;
t214 = (-pkin(2) * t237 - pkin(9) * t234 - pkin(1)) * t247;
t240 = -t233 * t212 + t253 * t214;
t197 = -t221 * pkin(3) - t216 * pkin(10) + t240;
t215 = -t253 * t226 + t233 * t243;
t249 = t253 * t212 + t233 * t214;
t200 = -t215 * pkin(10) + t249;
t232 = sin(qJ(4));
t236 = cos(qJ(4));
t250 = t232 * t197 + t236 * t200;
t190 = -t218 * qJ(5) + t250;
t205 = t236 * t215 + t232 * t216;
t206 = -t232 * t215 + t236 * t216;
t239 = -pkin(8) * t243 + t237 * t245;
t211 = -t226 * pkin(2) - t239;
t207 = t215 * pkin(3) + t211;
t193 = t205 * pkin(4) - t206 * qJ(5) + t207;
t228 = sin(pkin(12));
t186 = t252 * t190 + t228 * t193;
t244 = t237 * t251;
t185 = -t228 * t190 + t252 * t193;
t241 = t236 * t197 - t232 * t200;
t189 = t218 * pkin(4) + qJD(5) - t241;
t235 = cos(qJ(6));
t231 = sin(qJ(6));
t204 = qJD(6) + t205;
t203 = t252 * t206 - t228 * t218;
t202 = t228 * t206 + t252 * t218;
t195 = -t231 * t202 + t235 * t203;
t194 = t235 * t202 + t231 * t203;
t187 = t202 * pkin(5) + t189;
t184 = -t202 * pkin(11) + t186;
t183 = t205 * pkin(5) - t203 * pkin(11) + t185;
t1 = [t238 / 0.2e1, 0, 0, t234 ^ 2 * t251 / 0.2e1, t234 * t244, t226 * t243, t226 * t242, t226 ^ 2 / 0.2e1, pkin(1) * t244 + t239 * t226, -pkin(1) * t234 * t251 - t248 * t226, t216 ^ 2 / 0.2e1, -t216 * t215, -t216 * t221, t215 * t221, t221 ^ 2 / 0.2e1, t211 * t215 - t221 * t240, t211 * t216 + t249 * t221, t206 ^ 2 / 0.2e1, -t206 * t205, -t206 * t218, t205 * t218, t218 ^ 2 / 0.2e1, t207 * t205 - t218 * t241, t207 * t206 + t250 * t218, t185 * t205 + t189 * t202, -t186 * t205 + t189 * t203, -t185 * t203 - t186 * t202, t186 ^ 2 / 0.2e1 + t185 ^ 2 / 0.2e1 + t189 ^ 2 / 0.2e1, t195 ^ 2 / 0.2e1, -t195 * t194, t195 * t204, -t194 * t204, t204 ^ 2 / 0.2e1 (t235 * t183 - t231 * t184) * t204 + t187 * t194 -(t231 * t183 + t235 * t184) * t204 + t187 * t195;];
T_reg  = t1;
