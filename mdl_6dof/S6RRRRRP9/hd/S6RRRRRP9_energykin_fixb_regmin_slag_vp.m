% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRP9
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
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRP9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:11:38
% EndTime: 2019-03-10 02:11:38
% DurationCPUTime: 0.16s
% Computational Cost: add. (657->52), mult. (1503->108), div. (0->0), fcn. (1189->10), ass. (0->48)
t231 = cos(qJ(5));
t206 = sin(pkin(6));
t215 = qJD(1) ^ 2;
t230 = t206 ^ 2 * t215;
t224 = cos(pkin(6)) * qJD(1);
t204 = qJD(2) + t224;
t210 = sin(qJ(3));
t213 = cos(qJ(3));
t211 = sin(qJ(2));
t225 = qJD(1) * t206;
t221 = t211 * t225;
t196 = t210 * t204 + t213 * t221;
t214 = cos(qJ(2));
t220 = t214 * t225;
t199 = -qJD(3) + t220;
t209 = sin(qJ(4));
t212 = cos(qJ(4));
t186 = t212 * t196 - t209 * t199;
t195 = -t213 * t204 + t210 * t221;
t194 = qJD(4) + t195;
t223 = pkin(1) * t224;
t216 = -pkin(8) * t221 + t214 * t223;
t190 = -t204 * pkin(2) - t216;
t180 = pkin(3) * t195 - pkin(10) * t196 + t190;
t226 = pkin(8) * t220 + t211 * t223;
t191 = t204 * pkin(9) + t226;
t193 = (-pkin(2) * t214 - pkin(9) * t211 - pkin(1)) * t225;
t227 = t213 * t191 + t210 * t193;
t183 = -pkin(10) * t199 + t227;
t218 = t212 * t180 - t183 * t209;
t171 = pkin(4) * t194 - pkin(11) * t186 + t218;
t185 = t196 * t209 + t212 * t199;
t228 = t209 * t180 + t212 * t183;
t174 = -pkin(11) * t185 + t228;
t208 = sin(qJ(5));
t229 = t208 * t171 + t231 * t174;
t222 = t214 * t230;
t219 = t231 * t171 - t208 * t174;
t217 = -t210 * t191 + t193 * t213;
t182 = pkin(3) * t199 - t217;
t175 = pkin(4) * t185 + t182;
t192 = qJD(5) + t194;
t177 = -t208 * t185 + t231 * t186;
t176 = t231 * t185 + t208 * t186;
t172 = pkin(5) * t176 + qJD(6) + t175;
t168 = -qJ(6) * t176 + t229;
t167 = pkin(5) * t192 - qJ(6) * t177 + t219;
t1 = [t215 / 0.2e1, 0, 0, t211 ^ 2 * t230 / 0.2e1, t211 * t222, t204 * t221, t204 * t220, t204 ^ 2 / 0.2e1, pkin(1) * t222 + t216 * t204, -pkin(1) * t211 * t230 - t226 * t204, t196 ^ 2 / 0.2e1, -t196 * t195, -t196 * t199, t195 * t199, t199 ^ 2 / 0.2e1, t190 * t195 - t217 * t199, t190 * t196 + t227 * t199, t186 ^ 2 / 0.2e1, -t186 * t185, t186 * t194, -t185 * t194, t194 ^ 2 / 0.2e1, t182 * t185 + t218 * t194, t182 * t186 - t228 * t194, t177 ^ 2 / 0.2e1, -t177 * t176, t177 * t192, -t176 * t192, t192 ^ 2 / 0.2e1, t175 * t176 + t219 * t192, t175 * t177 - t229 * t192, -t167 * t177 - t168 * t176, t168 ^ 2 / 0.2e1 + t167 ^ 2 / 0.2e1 + t172 ^ 2 / 0.2e1;];
T_reg  = t1;
