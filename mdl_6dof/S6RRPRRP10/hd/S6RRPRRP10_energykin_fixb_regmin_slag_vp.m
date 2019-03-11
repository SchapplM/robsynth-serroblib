% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:42:31
% EndTime: 2019-03-09 12:42:31
% DurationCPUTime: 0.17s
% Computational Cost: add. (831->55), mult. (1984->113), div. (0->0), fcn. (1563->10), ass. (0->47)
t218 = sin(pkin(6));
t227 = qJD(1) ^ 2;
t240 = t218 ^ 2 * t227;
t226 = cos(qJ(2));
t236 = qJD(1) * t218;
t231 = t226 * t236;
t210 = -qJD(4) + t231;
t235 = cos(pkin(6)) * qJD(1);
t215 = qJD(2) + t235;
t223 = sin(qJ(2));
t234 = pkin(1) * t235;
t237 = pkin(8) * t231 + t223 * t234;
t203 = t215 * qJ(3) + t237;
t205 = (-pkin(2) * t226 - qJ(3) * t223 - pkin(1)) * t236;
t217 = sin(pkin(11));
t219 = cos(pkin(11));
t193 = -t217 * t203 + t219 * t205;
t232 = t223 * t236;
t207 = t217 * t215 + t219 * t232;
t187 = -pkin(3) * t231 - t207 * pkin(9) + t193;
t194 = t219 * t203 + t217 * t205;
t206 = -t219 * t215 + t217 * t232;
t190 = -t206 * pkin(9) + t194;
t222 = sin(qJ(4));
t225 = cos(qJ(4));
t238 = t222 * t187 + t225 * t190;
t183 = -t210 * pkin(10) + t238;
t228 = -pkin(8) * t232 + t226 * t234;
t200 = -t215 * pkin(2) + qJD(3) - t228;
t196 = t206 * pkin(3) + t200;
t197 = t225 * t206 + t222 * t207;
t198 = -t222 * t206 + t225 * t207;
t185 = t197 * pkin(4) - t198 * pkin(10) + t196;
t221 = sin(qJ(5));
t224 = cos(qJ(5));
t239 = t224 * t183 + t221 * t185;
t233 = t226 * t240;
t230 = t225 * t187 - t222 * t190;
t229 = -t221 * t183 + t224 * t185;
t182 = t210 * pkin(4) - t230;
t195 = qJD(5) + t197;
t192 = t224 * t198 - t221 * t210;
t191 = t221 * t198 + t224 * t210;
t180 = t191 * pkin(5) - t192 * qJ(6) + t182;
t179 = t195 * qJ(6) + t239;
t178 = -t195 * pkin(5) + qJD(6) - t229;
t1 = [t227 / 0.2e1, 0, 0, t223 ^ 2 * t240 / 0.2e1, t223 * t233, t215 * t232, t215 * t231, t215 ^ 2 / 0.2e1, pkin(1) * t233 + t228 * t215, -pkin(1) * t223 * t240 - t237 * t215, -t193 * t231 + t200 * t206, t194 * t231 + t200 * t207, -t193 * t207 - t194 * t206, t194 ^ 2 / 0.2e1 + t193 ^ 2 / 0.2e1 + t200 ^ 2 / 0.2e1, t198 ^ 2 / 0.2e1, -t198 * t197, -t198 * t210, t197 * t210, t210 ^ 2 / 0.2e1, t196 * t197 - t230 * t210, t196 * t198 + t238 * t210, t192 ^ 2 / 0.2e1, -t192 * t191, t192 * t195, -t191 * t195, t195 ^ 2 / 0.2e1, t182 * t191 + t229 * t195, t182 * t192 - t239 * t195, -t178 * t195 + t180 * t191, t178 * t192 - t179 * t191, t179 * t195 - t180 * t192, t179 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1 + t178 ^ 2 / 0.2e1;];
T_reg  = t1;
