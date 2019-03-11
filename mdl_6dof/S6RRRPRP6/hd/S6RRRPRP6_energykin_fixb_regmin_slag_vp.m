% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRP6
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
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:01:08
% EndTime: 2019-03-09 17:01:08
% DurationCPUTime: 0.15s
% Computational Cost: add. (632->51), mult. (1506->105), div. (0->0), fcn. (1175->10), ass. (0->47)
t238 = cos(qJ(3));
t237 = cos(qJ(5));
t216 = sin(pkin(6));
t223 = qJD(1) ^ 2;
t236 = t216 ^ 2 * t223;
t231 = cos(pkin(6)) * qJD(1);
t213 = qJD(2) + t231;
t220 = sin(qJ(3));
t221 = sin(qJ(2));
t232 = qJD(1) * t216;
t228 = t221 * t232;
t205 = t220 * t213 + t238 * t228;
t222 = cos(qJ(2));
t227 = t222 * t232;
t208 = -qJD(3) + t227;
t230 = pkin(1) * t231;
t233 = pkin(8) * t227 + t221 * t230;
t201 = t213 * pkin(9) + t233;
t203 = (-pkin(2) * t222 - pkin(9) * t221 - pkin(1)) * t232;
t225 = -t220 * t201 + t238 * t203;
t187 = -t208 * pkin(3) - t205 * qJ(4) + t225;
t204 = -t238 * t213 + t220 * t228;
t234 = t238 * t201 + t220 * t203;
t190 = -t204 * qJ(4) + t234;
t215 = sin(pkin(11));
t217 = cos(pkin(11));
t182 = t215 * t187 + t217 * t190;
t180 = -t208 * pkin(10) + t182;
t224 = -pkin(8) * t228 + t222 * t230;
t200 = -t213 * pkin(2) - t224;
t194 = t204 * pkin(3) + qJD(4) + t200;
t195 = -t217 * t204 - t215 * t205;
t196 = -t215 * t204 + t217 * t205;
t185 = -t195 * pkin(4) - t196 * pkin(10) + t194;
t219 = sin(qJ(5));
t235 = t237 * t180 + t219 * t185;
t229 = t222 * t236;
t226 = -t219 * t180 + t237 * t185;
t181 = t217 * t187 - t215 * t190;
t179 = t208 * pkin(4) - t181;
t193 = qJD(5) - t195;
t192 = t237 * t196 - t219 * t208;
t191 = t219 * t196 + t237 * t208;
t177 = t191 * pkin(5) + qJD(6) + t179;
t176 = -t191 * qJ(6) + t235;
t175 = t193 * pkin(5) - t192 * qJ(6) + t226;
t1 = [t223 / 0.2e1, 0, 0, t221 ^ 2 * t236 / 0.2e1, t221 * t229, t213 * t228, t213 * t227, t213 ^ 2 / 0.2e1, pkin(1) * t229 + t224 * t213, -pkin(1) * t221 * t236 - t233 * t213, t205 ^ 2 / 0.2e1, -t205 * t204, -t205 * t208, t204 * t208, t208 ^ 2 / 0.2e1, t200 * t204 - t225 * t208, t200 * t205 + t234 * t208, -t181 * t196 + t182 * t195, t182 ^ 2 / 0.2e1 + t181 ^ 2 / 0.2e1 + t194 ^ 2 / 0.2e1, t192 ^ 2 / 0.2e1, -t192 * t191, t192 * t193, -t191 * t193, t193 ^ 2 / 0.2e1, t179 * t191 + t226 * t193, t179 * t192 - t235 * t193, -t175 * t192 - t176 * t191, t176 ^ 2 / 0.2e1 + t175 ^ 2 / 0.2e1 + t177 ^ 2 / 0.2e1;];
T_reg  = t1;
