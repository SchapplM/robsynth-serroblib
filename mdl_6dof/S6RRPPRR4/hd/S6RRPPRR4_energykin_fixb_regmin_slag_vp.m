% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:02
% EndTime: 2019-03-09 09:05:02
% DurationCPUTime: 0.12s
% Computational Cost: add. (439->55), mult. (1243->110), div. (0->0), fcn. (953->10), ass. (0->48)
t237 = pkin(3) + pkin(9);
t215 = sin(pkin(6));
t224 = qJD(1) ^ 2;
t236 = t215 ^ 2 * t224;
t214 = sin(pkin(11));
t216 = cos(pkin(11));
t220 = sin(qJ(2));
t223 = cos(qJ(2));
t233 = qJD(1) * t215;
t205 = (t214 * t223 + t216 * t220) * t233;
t232 = cos(pkin(6)) * qJD(1);
t212 = qJD(2) + t232;
t231 = pkin(1) * t232;
t211 = t223 * t231;
t229 = t220 * t233;
t197 = t212 * pkin(2) + t211 + (-pkin(8) - qJ(3)) * t229;
t228 = t223 * t233;
t234 = pkin(8) * t228 + t220 * t231;
t200 = qJ(3) * t228 + t234;
t188 = t216 * t197 - t214 * t200;
t227 = qJD(4) - t188;
t182 = t205 * pkin(4) - t237 * t212 + t227;
t204 = t214 * t229 - t216 * t228;
t206 = qJD(3) + (-pkin(2) * t223 - pkin(1)) * t233;
t225 = -t205 * qJ(4) + t206;
t185 = t237 * t204 + t225;
t219 = sin(qJ(5));
t222 = cos(qJ(5));
t235 = t219 * t182 + t222 * t185;
t189 = t214 * t197 + t216 * t200;
t230 = t223 * t236;
t187 = -t212 * qJ(4) - t189;
t195 = -t222 * t204 + t219 * t212;
t226 = t222 * t182 - t219 * t185;
t183 = -t204 * pkin(4) - t187;
t221 = cos(qJ(6));
t218 = sin(qJ(6));
t203 = qJD(5) + t205;
t196 = t219 * t204 + t222 * t212;
t194 = qJD(6) + t195;
t192 = t204 * pkin(3) + t225;
t191 = t221 * t196 + t218 * t203;
t190 = t218 * t196 - t221 * t203;
t186 = -t212 * pkin(3) + t227;
t180 = t195 * pkin(5) - t196 * pkin(10) + t183;
t179 = t203 * pkin(10) + t235;
t178 = -t203 * pkin(5) - t226;
t1 = [t224 / 0.2e1, 0, 0, t220 ^ 2 * t236 / 0.2e1, t220 * t230, t212 * t229, t212 * t228, t212 ^ 2 / 0.2e1, pkin(1) * t230 + (-pkin(8) * t229 + t211) * t212, -pkin(1) * t220 * t236 - t234 * t212, -t188 * t205 - t189 * t204, t189 ^ 2 / 0.2e1 + t188 ^ 2 / 0.2e1 + t206 ^ 2 / 0.2e1, t186 * t205 + t187 * t204, t186 * t212 - t192 * t204, -t187 * t212 - t192 * t205, t192 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1 + t186 ^ 2 / 0.2e1, t196 ^ 2 / 0.2e1, -t196 * t195, t196 * t203, -t195 * t203, t203 ^ 2 / 0.2e1, t183 * t195 + t226 * t203, t183 * t196 - t235 * t203, t191 ^ 2 / 0.2e1, -t191 * t190, t191 * t194, -t190 * t194, t194 ^ 2 / 0.2e1 (-t218 * t179 + t221 * t180) * t194 + t178 * t190 -(t221 * t179 + t218 * t180) * t194 + t178 * t191;];
T_reg  = t1;
