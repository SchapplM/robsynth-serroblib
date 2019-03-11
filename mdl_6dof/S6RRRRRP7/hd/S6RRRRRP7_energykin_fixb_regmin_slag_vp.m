% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP7
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
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:43:08
% EndTime: 2019-03-10 01:43:09
% DurationCPUTime: 0.16s
% Computational Cost: add. (646->52), mult. (1487->108), div. (0->0), fcn. (1179->10), ass. (0->48)
t237 = cos(qJ(3));
t236 = cos(qJ(5));
t212 = sin(pkin(6));
t220 = qJD(1) ^ 2;
t235 = t212 ^ 2 * t220;
t219 = cos(qJ(2));
t230 = qJD(1) * t212;
t225 = t219 * t230;
t205 = -qJD(3) + t225;
t202 = -qJD(4) + t205;
t229 = cos(pkin(6)) * qJD(1);
t210 = qJD(2) + t229;
t216 = sin(qJ(3));
t217 = sin(qJ(2));
t226 = t217 * t230;
t200 = t216 * t210 + t237 * t226;
t228 = pkin(1) * t229;
t231 = pkin(8) * t225 + t217 * t228;
t196 = t210 * pkin(9) + t231;
t198 = (-pkin(2) * t219 - pkin(9) * t217 - pkin(1)) * t230;
t222 = -t216 * t196 + t237 * t198;
t182 = -t205 * pkin(3) - t200 * pkin(10) + t222;
t199 = -t237 * t210 + t216 * t226;
t232 = t237 * t196 + t216 * t198;
t185 = -t199 * pkin(10) + t232;
t215 = sin(qJ(4));
t218 = cos(qJ(4));
t233 = t215 * t182 + t218 * t185;
t177 = -t202 * pkin(11) + t233;
t189 = t218 * t199 + t215 * t200;
t190 = -t215 * t199 + t218 * t200;
t221 = -pkin(8) * t226 + t219 * t228;
t195 = -t210 * pkin(2) - t221;
t191 = t199 * pkin(3) + t195;
t180 = t189 * pkin(4) - t190 * pkin(11) + t191;
t214 = sin(qJ(5));
t234 = t236 * t177 + t214 * t180;
t227 = t219 * t235;
t224 = -t214 * t177 + t236 * t180;
t223 = t218 * t182 - t215 * t185;
t176 = t202 * pkin(4) - t223;
t188 = qJD(5) + t189;
t187 = t236 * t190 - t214 * t202;
t186 = t214 * t190 + t236 * t202;
t174 = t186 * pkin(5) + qJD(6) + t176;
t173 = -t186 * qJ(6) + t234;
t172 = t188 * pkin(5) - t187 * qJ(6) + t224;
t1 = [t220 / 0.2e1, 0, 0, t217 ^ 2 * t235 / 0.2e1, t217 * t227, t210 * t226, t210 * t225, t210 ^ 2 / 0.2e1, pkin(1) * t227 + t221 * t210, -pkin(1) * t217 * t235 - t231 * t210, t200 ^ 2 / 0.2e1, -t200 * t199, -t200 * t205, t199 * t205, t205 ^ 2 / 0.2e1, t195 * t199 - t222 * t205, t195 * t200 + t232 * t205, t190 ^ 2 / 0.2e1, -t190 * t189, -t190 * t202, t189 * t202, t202 ^ 2 / 0.2e1, t191 * t189 - t223 * t202, t191 * t190 + t233 * t202, t187 ^ 2 / 0.2e1, -t187 * t186, t187 * t188, -t186 * t188, t188 ^ 2 / 0.2e1, t176 * t186 + t224 * t188, t176 * t187 - t234 * t188, -t172 * t187 - t173 * t186, t173 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1 + t174 ^ 2 / 0.2e1;];
T_reg  = t1;
