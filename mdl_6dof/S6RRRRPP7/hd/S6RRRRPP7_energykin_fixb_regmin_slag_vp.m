% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:26:30
% EndTime: 2019-03-09 21:26:30
% DurationCPUTime: 0.18s
% Computational Cost: add. (866->53), mult. (1961->109), div. (0->0), fcn. (1525->10), ass. (0->47)
t247 = cos(qJ(4));
t225 = sin(pkin(6));
t233 = qJD(1) ^ 2;
t246 = t225 ^ 2 * t233;
t241 = cos(pkin(6)) * qJD(1);
t222 = qJD(2) + t241;
t229 = sin(qJ(3));
t231 = cos(qJ(3));
t230 = sin(qJ(2));
t242 = qJD(1) * t225;
t238 = t230 * t242;
t214 = t222 * t229 + t231 * t238;
t232 = cos(qJ(2));
t237 = t232 * t242;
t217 = -qJD(3) + t237;
t228 = sin(qJ(4));
t205 = t247 * t214 - t228 * t217;
t213 = -t231 * t222 + t229 * t238;
t212 = qJD(4) + t213;
t240 = pkin(1) * t241;
t234 = -pkin(8) * t238 + t232 * t240;
t209 = -t222 * pkin(2) - t234;
t199 = t213 * pkin(3) - t214 * pkin(10) + t209;
t243 = pkin(8) * t237 + t230 * t240;
t210 = pkin(9) * t222 + t243;
t211 = (-pkin(2) * t232 - pkin(9) * t230 - pkin(1)) * t242;
t244 = t231 * t210 + t229 * t211;
t202 = -pkin(10) * t217 + t244;
t236 = t247 * t199 - t228 * t202;
t191 = t212 * pkin(4) - t205 * qJ(5) + t236;
t204 = t214 * t228 + t247 * t217;
t245 = t228 * t199 + t247 * t202;
t193 = -qJ(5) * t204 + t245;
t224 = sin(pkin(11));
t226 = cos(pkin(11));
t188 = t224 * t191 + t226 * t193;
t239 = t232 * t246;
t235 = -t229 * t210 + t211 * t231;
t187 = t191 * t226 - t193 * t224;
t201 = t217 * pkin(3) - t235;
t194 = t204 * pkin(4) + qJD(5) + t201;
t196 = -t204 * t224 + t205 * t226;
t195 = t226 * t204 + t205 * t224;
t189 = t195 * pkin(5) - t196 * qJ(6) + t194;
t186 = qJ(6) * t212 + t188;
t185 = -pkin(5) * t212 + qJD(6) - t187;
t1 = [t233 / 0.2e1, 0, 0, t230 ^ 2 * t246 / 0.2e1, t230 * t239, t222 * t238, t222 * t237, t222 ^ 2 / 0.2e1, pkin(1) * t239 + t234 * t222, -pkin(1) * t230 * t246 - t243 * t222, t214 ^ 2 / 0.2e1, -t214 * t213, -t214 * t217, t213 * t217, t217 ^ 2 / 0.2e1, t209 * t213 - t235 * t217, t209 * t214 + t244 * t217, t205 ^ 2 / 0.2e1, -t205 * t204, t205 * t212, -t204 * t212, t212 ^ 2 / 0.2e1, t201 * t204 + t236 * t212, t201 * t205 - t245 * t212, -t187 * t196 - t188 * t195, t188 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1 + t194 ^ 2 / 0.2e1, -t185 * t212 + t189 * t195, t185 * t196 - t186 * t195, t186 * t212 - t189 * t196, t186 ^ 2 / 0.2e1 + t189 ^ 2 / 0.2e1 + t185 ^ 2 / 0.2e1;];
T_reg  = t1;
