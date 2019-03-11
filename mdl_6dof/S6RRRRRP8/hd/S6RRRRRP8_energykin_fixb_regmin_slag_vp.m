% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP8
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:55:48
% EndTime: 2019-03-10 01:55:48
% DurationCPUTime: 0.14s
% Computational Cost: add. (849->54), mult. (1924->112), div. (0->0), fcn. (1531->10), ass. (0->48)
t244 = cos(qJ(3));
t219 = sin(pkin(6));
t228 = qJD(1) ^ 2;
t243 = t219 ^ 2 * t228;
t227 = cos(qJ(2));
t238 = qJD(1) * t219;
t233 = t227 * t238;
t212 = -qJD(3) + t233;
t209 = -qJD(4) + t212;
t237 = cos(pkin(6)) * qJD(1);
t217 = qJD(2) + t237;
t223 = sin(qJ(3));
t224 = sin(qJ(2));
t234 = t224 * t238;
t207 = t223 * t217 + t244 * t234;
t236 = pkin(1) * t237;
t239 = pkin(8) * t233 + t224 * t236;
t203 = t217 * pkin(9) + t239;
t205 = (-pkin(2) * t227 - pkin(9) * t224 - pkin(1)) * t238;
t231 = -t223 * t203 + t244 * t205;
t189 = -t212 * pkin(3) - t207 * pkin(10) + t231;
t206 = -t244 * t217 + t223 * t234;
t240 = t244 * t203 + t223 * t205;
t192 = -t206 * pkin(10) + t240;
t222 = sin(qJ(4));
t226 = cos(qJ(4));
t241 = t222 * t189 + t226 * t192;
t185 = -t209 * pkin(11) + t241;
t196 = t226 * t206 + t222 * t207;
t197 = -t222 * t206 + t226 * t207;
t229 = -pkin(8) * t234 + t227 * t236;
t202 = -t217 * pkin(2) - t229;
t198 = t206 * pkin(3) + t202;
t187 = t196 * pkin(4) - t197 * pkin(11) + t198;
t221 = sin(qJ(5));
t225 = cos(qJ(5));
t242 = t225 * t185 + t221 * t187;
t235 = t227 * t243;
t232 = t226 * t189 - t222 * t192;
t230 = -t221 * t185 + t225 * t187;
t184 = t209 * pkin(4) - t232;
t195 = qJD(5) + t196;
t194 = t225 * t197 - t221 * t209;
t193 = t221 * t197 + t225 * t209;
t182 = t193 * pkin(5) - t194 * qJ(6) + t184;
t181 = t195 * qJ(6) + t242;
t180 = -t195 * pkin(5) + qJD(6) - t230;
t1 = [t228 / 0.2e1, 0, 0, t224 ^ 2 * t243 / 0.2e1, t224 * t235, t217 * t234, t217 * t233, t217 ^ 2 / 0.2e1, pkin(1) * t235 + t229 * t217, -pkin(1) * t224 * t243 - t239 * t217, t207 ^ 2 / 0.2e1, -t207 * t206, -t207 * t212, t206 * t212, t212 ^ 2 / 0.2e1, t202 * t206 - t231 * t212, t202 * t207 + t240 * t212, t197 ^ 2 / 0.2e1, -t197 * t196, -t197 * t209, t196 * t209, t209 ^ 2 / 0.2e1, t198 * t196 - t232 * t209, t198 * t197 + t241 * t209, t194 ^ 2 / 0.2e1, -t194 * t193, t194 * t195, -t193 * t195, t195 ^ 2 / 0.2e1, t184 * t193 + t230 * t195, t184 * t194 - t242 * t195, -t180 * t195 + t182 * t193, t180 * t194 - t181 * t193, t181 * t195 - t182 * t194, t181 ^ 2 / 0.2e1 + t182 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1;];
T_reg  = t1;
