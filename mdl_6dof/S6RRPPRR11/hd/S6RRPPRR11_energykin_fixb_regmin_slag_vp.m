% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPPRR11_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:42:43
% EndTime: 2019-03-09 09:42:43
% DurationCPUTime: 0.17s
% Computational Cost: add. (506->58), mult. (1207->119), div. (0->0), fcn. (890->10), ass. (0->49)
t228 = -pkin(2) - qJ(4);
t206 = sin(pkin(6));
t215 = qJD(1) ^ 2;
t227 = t206 ^ 2 * t215;
t211 = sin(qJ(2));
t224 = qJD(1) * t206;
t219 = t211 * t224;
t199 = pkin(8) * t219;
t208 = cos(pkin(6));
t223 = t208 * qJD(1);
t203 = qJD(2) + t223;
t214 = cos(qJ(2));
t186 = qJD(3) + t199 + t228 * t203 + (-pkin(1) * t208 * t214 + pkin(3) * t206 * t211) * qJD(1);
t218 = -qJ(3) * t211 - pkin(1);
t189 = (t228 * t214 + t218) * t224;
t205 = sin(pkin(11));
t207 = cos(pkin(11));
t176 = t207 * t186 - t205 * t189;
t220 = t214 * t224;
t195 = t207 * t203 - t205 * t220;
t173 = pkin(4) * t219 - t195 * pkin(9) + t176;
t177 = t205 * t186 + t207 * t189;
t194 = t205 * t203 + t207 * t220;
t175 = -pkin(9) * t194 + t177;
t210 = sin(qJ(5));
t213 = cos(qJ(5));
t226 = t210 * t173 + t213 * t175;
t222 = pkin(1) * t223;
t225 = pkin(8) * t220 + t211 * t222;
t221 = t214 * t227;
t191 = -t203 * qJ(3) - t225;
t182 = t213 * t194 + t210 * t195;
t187 = pkin(3) * t220 + qJD(4) - t191;
t217 = t213 * t173 - t210 * t175;
t216 = t214 * t222 - t199;
t180 = t194 * pkin(4) + t187;
t212 = cos(qJ(6));
t209 = sin(qJ(6));
t197 = qJD(5) + t219;
t193 = (-pkin(2) * t214 + t218) * t224;
t190 = -t203 * pkin(2) + qJD(3) - t216;
t183 = -t210 * t194 + t213 * t195;
t181 = qJD(6) + t182;
t179 = t212 * t183 + t209 * t197;
t178 = t209 * t183 - t212 * t197;
t171 = pkin(5) * t182 - pkin(10) * t183 + t180;
t170 = pkin(10) * t197 + t226;
t169 = -t197 * pkin(5) - t217;
t1 = [t215 / 0.2e1, 0, 0, t211 ^ 2 * t227 / 0.2e1, t211 * t221, t203 * t219, t203 * t220, t203 ^ 2 / 0.2e1, pkin(1) * t221 + t216 * t203, -pkin(1) * t211 * t227 - t225 * t203 (t190 * t211 - t191 * t214) * t224, t190 * t203 + t193 * t220, -t191 * t203 - t193 * t219, t193 ^ 2 / 0.2e1 + t191 ^ 2 / 0.2e1 + t190 ^ 2 / 0.2e1, t176 * t219 + t187 * t194, -t177 * t219 + t187 * t195, -t176 * t195 - t177 * t194, t177 ^ 2 / 0.2e1 + t176 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1, t183 ^ 2 / 0.2e1, -t183 * t182, t183 * t197, -t182 * t197, t197 ^ 2 / 0.2e1, t180 * t182 + t217 * t197, t180 * t183 - t226 * t197, t179 ^ 2 / 0.2e1, -t179 * t178, t179 * t181, -t178 * t181, t181 ^ 2 / 0.2e1 (-t209 * t170 + t212 * t171) * t181 + t169 * t178 -(t212 * t170 + t209 * t171) * t181 + t169 * t179;];
T_reg  = t1;
