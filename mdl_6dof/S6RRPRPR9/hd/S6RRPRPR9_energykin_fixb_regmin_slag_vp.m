% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% T_reg [1x32]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR9_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:01:08
% EndTime: 2019-03-09 11:01:08
% DurationCPUTime: 0.20s
% Computational Cost: add. (943->59), mult. (2258->121), div. (0->0), fcn. (1797->12), ass. (0->51)
t249 = cos(pkin(12));
t228 = sin(pkin(6));
t237 = qJD(1) ^ 2;
t248 = t228 ^ 2 * t237;
t236 = cos(qJ(2));
t245 = qJD(1) * t228;
t240 = t236 * t245;
t219 = -qJD(4) + t240;
t244 = cos(pkin(6)) * qJD(1);
t224 = qJD(2) + t244;
t233 = sin(qJ(2));
t243 = pkin(1) * t244;
t246 = pkin(8) * t240 + t233 * t243;
t212 = t224 * qJ(3) + t246;
t214 = (-pkin(2) * t236 - qJ(3) * t233 - pkin(1)) * t245;
t227 = sin(pkin(11));
t229 = cos(pkin(11));
t202 = -t227 * t212 + t229 * t214;
t241 = t233 * t245;
t216 = t227 * t224 + t229 * t241;
t195 = -pkin(3) * t240 - t216 * pkin(9) + t202;
t203 = t229 * t212 + t227 * t214;
t215 = -t229 * t224 + t227 * t241;
t198 = -t215 * pkin(9) + t203;
t232 = sin(qJ(4));
t235 = cos(qJ(4));
t247 = t232 * t195 + t235 * t198;
t188 = -t219 * qJ(5) + t247;
t238 = -pkin(8) * t241 + t236 * t243;
t209 = -t224 * pkin(2) + qJD(3) - t238;
t205 = t215 * pkin(3) + t209;
t206 = t235 * t215 + t232 * t216;
t207 = -t232 * t215 + t235 * t216;
t191 = t206 * pkin(4) - t207 * qJ(5) + t205;
t226 = sin(pkin(12));
t184 = t249 * t188 + t226 * t191;
t242 = t236 * t248;
t183 = -t226 * t188 + t249 * t191;
t239 = t235 * t195 - t232 * t198;
t187 = t219 * pkin(4) + qJD(5) - t239;
t234 = cos(qJ(6));
t231 = sin(qJ(6));
t204 = qJD(6) + t206;
t201 = t249 * t207 - t226 * t219;
t200 = t226 * t207 + t249 * t219;
t193 = -t231 * t200 + t234 * t201;
t192 = t234 * t200 + t231 * t201;
t185 = t200 * pkin(5) + t187;
t182 = -t200 * pkin(10) + t184;
t181 = t206 * pkin(5) - t201 * pkin(10) + t183;
t1 = [t237 / 0.2e1, 0, 0, t233 ^ 2 * t248 / 0.2e1, t233 * t242, t224 * t241, t224 * t240, t224 ^ 2 / 0.2e1, pkin(1) * t242 + t238 * t224, -pkin(1) * t233 * t248 - t246 * t224, -t202 * t240 + t209 * t215, t203 * t240 + t209 * t216, -t202 * t216 - t203 * t215, t203 ^ 2 / 0.2e1 + t202 ^ 2 / 0.2e1 + t209 ^ 2 / 0.2e1, t207 ^ 2 / 0.2e1, -t207 * t206, -t207 * t219, t206 * t219, t219 ^ 2 / 0.2e1, t205 * t206 - t239 * t219, t205 * t207 + t247 * t219, t183 * t206 + t187 * t200, -t184 * t206 + t187 * t201, -t183 * t201 - t184 * t200, t184 ^ 2 / 0.2e1 + t183 ^ 2 / 0.2e1 + t187 ^ 2 / 0.2e1, t193 ^ 2 / 0.2e1, -t193 * t192, t193 * t204, -t192 * t204, t204 ^ 2 / 0.2e1 (t234 * t181 - t231 * t182) * t204 + t185 * t192 -(t231 * t181 + t234 * t182) * t204 + t185 * t193;];
T_reg  = t1;
