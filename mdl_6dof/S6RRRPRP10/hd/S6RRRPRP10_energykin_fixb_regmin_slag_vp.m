% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:36:55
% EndTime: 2019-03-09 17:36:55
% DurationCPUTime: 0.17s
% Computational Cost: add. (892->55), mult. (2044->113), div. (0->0), fcn. (1595->10), ass. (0->47)
t219 = sin(pkin(6));
t228 = qJD(1) ^ 2;
t241 = t219 ^ 2 * t228;
t236 = cos(pkin(6)) * qJD(1);
t216 = qJD(2) + t236;
t227 = cos(qJ(2));
t224 = sin(qJ(2));
t237 = qJD(1) * t219;
t233 = t224 * t237;
t235 = pkin(1) * t236;
t229 = -pkin(8) * t233 + t227 * t235;
t203 = -t216 * pkin(2) - t229;
t223 = sin(qJ(3));
t226 = cos(qJ(3));
t207 = -t226 * t216 + t223 * t233;
t208 = t223 * t216 + t226 * t233;
t193 = t207 * pkin(3) - t208 * qJ(4) + t203;
t232 = t227 * t237;
t211 = -qJD(3) + t232;
t238 = pkin(8) * t232 + t224 * t235;
t204 = t216 * pkin(9) + t238;
t205 = (-pkin(2) * t227 - pkin(9) * t224 - pkin(1)) * t237;
t239 = t226 * t204 + t223 * t205;
t196 = -t211 * qJ(4) + t239;
t218 = sin(pkin(11));
t220 = cos(pkin(11));
t186 = t220 * t193 - t218 * t196;
t199 = t220 * t208 - t218 * t211;
t183 = t207 * pkin(4) - t199 * pkin(10) + t186;
t187 = t218 * t193 + t220 * t196;
t198 = t218 * t208 + t220 * t211;
t185 = -t198 * pkin(10) + t187;
t222 = sin(qJ(5));
t225 = cos(qJ(5));
t240 = t222 * t183 + t225 * t185;
t234 = t227 * t241;
t231 = -t223 * t204 + t226 * t205;
t230 = t225 * t183 - t222 * t185;
t195 = t211 * pkin(3) + qJD(4) - t231;
t188 = t198 * pkin(4) + t195;
t206 = qJD(5) + t207;
t190 = -t222 * t198 + t225 * t199;
t189 = t225 * t198 + t222 * t199;
t181 = t189 * pkin(5) - t190 * qJ(6) + t188;
t180 = t206 * qJ(6) + t240;
t179 = -t206 * pkin(5) + qJD(6) - t230;
t1 = [t228 / 0.2e1, 0, 0, t224 ^ 2 * t241 / 0.2e1, t224 * t234, t216 * t233, t216 * t232, t216 ^ 2 / 0.2e1, pkin(1) * t234 + t229 * t216, -pkin(1) * t224 * t241 - t238 * t216, t208 ^ 2 / 0.2e1, -t208 * t207, -t208 * t211, t207 * t211, t211 ^ 2 / 0.2e1, t203 * t207 - t231 * t211, t203 * t208 + t239 * t211, t186 * t207 + t195 * t198, -t187 * t207 + t195 * t199, -t186 * t199 - t187 * t198, t187 ^ 2 / 0.2e1 + t186 ^ 2 / 0.2e1 + t195 ^ 2 / 0.2e1, t190 ^ 2 / 0.2e1, -t190 * t189, t190 * t206, -t189 * t206, t206 ^ 2 / 0.2e1, t188 * t189 + t230 * t206, t188 * t190 - t240 * t206, -t179 * t206 + t181 * t189, t179 * t190 - t180 * t189, t180 * t206 - t181 * t190, t180 ^ 2 / 0.2e1 + t181 ^ 2 / 0.2e1 + t179 ^ 2 / 0.2e1;];
T_reg  = t1;
