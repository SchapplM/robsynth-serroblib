% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR13_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:04:03
% EndTime: 2019-03-10 00:04:03
% DurationCPUTime: 0.13s
% Computational Cost: add. (612->55), mult. (1390->112), div. (0->0), fcn. (1082->10), ass. (0->49)
t243 = -pkin(4) - pkin(5);
t242 = cos(qJ(3));
t217 = sin(pkin(6));
t226 = qJD(1) ^ 2;
t241 = t217 ^ 2 * t226;
t235 = cos(pkin(6)) * qJD(1);
t215 = qJD(2) + t235;
t225 = cos(qJ(2));
t222 = sin(qJ(2));
t236 = qJD(1) * t217;
t232 = t222 * t236;
t234 = pkin(1) * t235;
t228 = -pkin(8) * t232 + t225 * t234;
t199 = -t215 * pkin(2) - t228;
t221 = sin(qJ(3));
t206 = -t242 * t215 + t221 * t232;
t207 = t221 * t215 + t242 * t232;
t187 = t206 * pkin(3) - t207 * pkin(10) + t199;
t231 = t225 * t236;
t210 = -qJD(3) + t231;
t237 = pkin(8) * t231 + t222 * t234;
t200 = t215 * pkin(9) + t237;
t204 = (-pkin(2) * t225 - pkin(9) * t222 - pkin(1)) * t236;
t238 = t242 * t200 + t221 * t204;
t191 = -t210 * pkin(10) + t238;
t220 = sin(qJ(4));
t224 = cos(qJ(4));
t240 = t220 * t187 + t224 * t191;
t239 = -t221 * t200 + t242 * t204;
t233 = t225 * t241;
t205 = qJD(4) + t206;
t182 = t205 * qJ(5) + t240;
t230 = t224 * t187 - t220 * t191;
t190 = t210 * pkin(3) - t239;
t229 = qJD(5) - t230;
t194 = t224 * t207 - t220 * t210;
t227 = t194 * qJ(5) - t190;
t223 = cos(qJ(6));
t219 = sin(qJ(6));
t202 = -qJD(6) + t205;
t193 = t220 * t207 + t224 * t210;
t185 = t219 * t193 + t223 * t194;
t184 = -t223 * t193 + t219 * t194;
t183 = t193 * pkin(4) - t227;
t181 = -t205 * pkin(4) + t229;
t180 = t243 * t193 + t227;
t179 = t193 * pkin(11) + t182;
t178 = -t194 * pkin(11) + t243 * t205 + t229;
t1 = [t226 / 0.2e1, 0, 0, t222 ^ 2 * t241 / 0.2e1, t222 * t233, t215 * t232, t215 * t231, t215 ^ 2 / 0.2e1, pkin(1) * t233 + t228 * t215, -pkin(1) * t222 * t241 - t237 * t215, t207 ^ 2 / 0.2e1, -t207 * t206, -t207 * t210, t206 * t210, t210 ^ 2 / 0.2e1, t199 * t206 - t239 * t210, t199 * t207 + t238 * t210, t194 ^ 2 / 0.2e1, -t194 * t193, t194 * t205, -t193 * t205, t205 ^ 2 / 0.2e1, t190 * t193 + t230 * t205, t190 * t194 - t240 * t205, -t181 * t205 + t183 * t193, t181 * t194 - t182 * t193, t182 * t205 - t183 * t194, t182 ^ 2 / 0.2e1 + t183 ^ 2 / 0.2e1 + t181 ^ 2 / 0.2e1, t185 ^ 2 / 0.2e1, -t185 * t184, -t185 * t202, t184 * t202, t202 ^ 2 / 0.2e1 -(t223 * t178 - t219 * t179) * t202 + t180 * t184 (t219 * t178 + t223 * t179) * t202 + t180 * t185;];
T_reg  = t1;
