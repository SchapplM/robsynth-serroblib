% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPPR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:17:54
% EndTime: 2019-03-09 16:17:54
% DurationCPUTime: 0.17s
% Computational Cost: add. (661->56), mult. (1542->113), div. (0->0), fcn. (1178->10), ass. (0->48)
t244 = -pkin(4) - pkin(5);
t243 = cos(qJ(3));
t221 = sin(pkin(6));
t229 = qJD(1) ^ 2;
t242 = t221 ^ 2 * t229;
t237 = cos(pkin(6)) * qJD(1);
t218 = qJD(2) + t237;
t228 = cos(qJ(2));
t226 = sin(qJ(2));
t238 = qJD(1) * t221;
t234 = t226 * t238;
t236 = pkin(1) * t237;
t231 = -pkin(8) * t234 + t228 * t236;
t203 = -t218 * pkin(2) - t231;
t225 = sin(qJ(3));
t209 = -t243 * t218 + t225 * t234;
t210 = t225 * t218 + t243 * t234;
t191 = t209 * pkin(3) - t210 * qJ(4) + t203;
t233 = t228 * t238;
t213 = -qJD(3) + t233;
t239 = pkin(8) * t233 + t226 * t236;
t204 = t218 * pkin(9) + t239;
t205 = (-pkin(2) * t228 - pkin(9) * t226 - pkin(1)) * t238;
t240 = t243 * t204 + t225 * t205;
t195 = -t213 * qJ(4) + t240;
t220 = sin(pkin(11));
t222 = cos(pkin(11));
t187 = t220 * t191 + t222 * t195;
t241 = -t225 * t204 + t243 * t205;
t235 = t228 * t242;
t184 = t209 * qJ(5) + t187;
t186 = t222 * t191 - t220 * t195;
t232 = qJD(5) - t186;
t194 = t213 * pkin(3) + qJD(4) - t241;
t198 = t222 * t210 - t220 * t213;
t230 = t198 * qJ(5) - t194;
t227 = cos(qJ(6));
t224 = sin(qJ(6));
t207 = -qJD(6) + t209;
t197 = t220 * t210 + t222 * t213;
t189 = t224 * t197 + t227 * t198;
t188 = -t227 * t197 + t224 * t198;
t185 = t197 * pkin(4) - t230;
t183 = -t209 * pkin(4) + t232;
t182 = t244 * t197 + t230;
t181 = t197 * pkin(10) + t184;
t180 = -t198 * pkin(10) + t244 * t209 + t232;
t1 = [t229 / 0.2e1, 0, 0, t226 ^ 2 * t242 / 0.2e1, t226 * t235, t218 * t234, t218 * t233, t218 ^ 2 / 0.2e1, pkin(1) * t235 + t231 * t218, -pkin(1) * t226 * t242 - t239 * t218, t210 ^ 2 / 0.2e1, -t210 * t209, -t210 * t213, t209 * t213, t213 ^ 2 / 0.2e1, t203 * t209 - t241 * t213, t203 * t210 + t240 * t213, t186 * t209 + t194 * t197, -t187 * t209 + t194 * t198, -t186 * t198 - t187 * t197, t187 ^ 2 / 0.2e1 + t186 ^ 2 / 0.2e1 + t194 ^ 2 / 0.2e1, -t183 * t209 + t185 * t197, t183 * t198 - t184 * t197, t184 * t209 - t185 * t198, t184 ^ 2 / 0.2e1 + t185 ^ 2 / 0.2e1 + t183 ^ 2 / 0.2e1, t189 ^ 2 / 0.2e1, -t189 * t188, -t189 * t207, t188 * t207, t207 ^ 2 / 0.2e1 -(t227 * t180 - t224 * t181) * t207 + t182 * t188 (t224 * t180 + t227 * t181) * t207 + t182 * t189;];
T_reg  = t1;
