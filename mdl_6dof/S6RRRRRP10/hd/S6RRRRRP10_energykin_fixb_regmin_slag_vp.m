% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:27:25
% EndTime: 2019-03-10 02:27:25
% DurationCPUTime: 0.14s
% Computational Cost: add. (845->54), mult. (1892->112), div. (0->0), fcn. (1499->10), ass. (0->48)
t243 = cos(qJ(4));
t218 = sin(pkin(6));
t227 = qJD(1) ^ 2;
t242 = t218 ^ 2 * t227;
t236 = cos(pkin(6)) * qJD(1);
t216 = qJD(2) + t236;
t222 = sin(qJ(3));
t225 = cos(qJ(3));
t223 = sin(qJ(2));
t237 = qJD(1) * t218;
t233 = t223 * t237;
t208 = t222 * t216 + t225 * t233;
t226 = cos(qJ(2));
t232 = t226 * t237;
t211 = -qJD(3) + t232;
t221 = sin(qJ(4));
t198 = t243 * t208 - t221 * t211;
t207 = -t225 * t216 + t222 * t233;
t206 = qJD(4) + t207;
t235 = pkin(1) * t236;
t228 = -pkin(8) * t233 + t226 * t235;
t202 = -t216 * pkin(2) - t228;
t192 = t207 * pkin(3) - t208 * pkin(10) + t202;
t238 = pkin(8) * t232 + t223 * t235;
t203 = t216 * pkin(9) + t238;
t205 = (-pkin(2) * t226 - pkin(9) * t223 - pkin(1)) * t237;
t239 = t225 * t203 + t222 * t205;
t195 = -t211 * pkin(10) + t239;
t231 = t243 * t192 - t221 * t195;
t184 = t206 * pkin(4) - t198 * pkin(11) + t231;
t197 = t221 * t208 + t243 * t211;
t240 = t221 * t192 + t243 * t195;
t186 = -t197 * pkin(11) + t240;
t220 = sin(qJ(5));
t224 = cos(qJ(5));
t241 = t220 * t184 + t224 * t186;
t234 = t226 * t242;
t230 = -t222 * t203 + t225 * t205;
t229 = t224 * t184 - t220 * t186;
t194 = t211 * pkin(3) - t230;
t187 = t197 * pkin(4) + t194;
t204 = qJD(5) + t206;
t189 = -t220 * t197 + t224 * t198;
t188 = t224 * t197 + t220 * t198;
t182 = t188 * pkin(5) - t189 * qJ(6) + t187;
t181 = t204 * qJ(6) + t241;
t180 = -t204 * pkin(5) + qJD(6) - t229;
t1 = [t227 / 0.2e1, 0, 0, t223 ^ 2 * t242 / 0.2e1, t223 * t234, t216 * t233, t216 * t232, t216 ^ 2 / 0.2e1, pkin(1) * t234 + t228 * t216, -pkin(1) * t223 * t242 - t238 * t216, t208 ^ 2 / 0.2e1, -t208 * t207, -t208 * t211, t207 * t211, t211 ^ 2 / 0.2e1, t202 * t207 - t230 * t211, t202 * t208 + t239 * t211, t198 ^ 2 / 0.2e1, -t198 * t197, t198 * t206, -t197 * t206, t206 ^ 2 / 0.2e1, t194 * t197 + t231 * t206, t194 * t198 - t240 * t206, t189 ^ 2 / 0.2e1, -t189 * t188, t189 * t204, -t188 * t204, t204 ^ 2 / 0.2e1, t187 * t188 + t229 * t204, t187 * t189 - t241 * t204, -t180 * t204 + t182 * t188, t180 * t189 - t181 * t188, t181 * t204 - t182 * t189, t181 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1;];
T_reg  = t1;
