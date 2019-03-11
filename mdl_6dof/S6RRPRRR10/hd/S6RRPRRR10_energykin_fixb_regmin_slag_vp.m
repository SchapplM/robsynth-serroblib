% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRR10_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:25:22
% EndTime: 2019-03-09 14:25:22
% DurationCPUTime: 0.17s
% Computational Cost: add. (793->58), mult. (1893->120), div. (0->0), fcn. (1537->12), ass. (0->52)
t251 = cos(qJ(5));
t250 = cos(pkin(12));
t227 = sin(pkin(6));
t236 = qJD(1) ^ 2;
t249 = t227 ^ 2 * t236;
t235 = cos(qJ(2));
t245 = qJD(1) * t227;
t240 = t235 * t245;
t219 = -qJD(4) + t240;
t244 = cos(pkin(6)) * qJD(1);
t224 = qJD(2) + t244;
t232 = sin(qJ(2));
t243 = pkin(1) * t244;
t246 = pkin(8) * t240 + t232 * t243;
t212 = t224 * qJ(3) + t246;
t214 = (-pkin(2) * t235 - qJ(3) * t232 - pkin(1)) * t245;
t226 = sin(pkin(12));
t201 = -t226 * t212 + t250 * t214;
t241 = t232 * t245;
t216 = t226 * t224 + t241 * t250;
t194 = -pkin(3) * t240 - t216 * pkin(9) + t201;
t202 = t250 * t212 + t226 * t214;
t215 = -t224 * t250 + t226 * t241;
t197 = -t215 * pkin(9) + t202;
t231 = sin(qJ(4));
t234 = cos(qJ(4));
t247 = t231 * t194 + t234 * t197;
t187 = -t219 * pkin(10) + t247;
t237 = -pkin(8) * t241 + t235 * t243;
t209 = -t224 * pkin(2) + qJD(3) - t237;
t205 = t215 * pkin(3) + t209;
t206 = t234 * t215 + t231 * t216;
t207 = -t231 * t215 + t234 * t216;
t190 = t206 * pkin(4) - t207 * pkin(10) + t205;
t230 = sin(qJ(5));
t248 = t187 * t251 + t230 * t190;
t242 = t235 * t249;
t239 = -t230 * t187 + t190 * t251;
t238 = t234 * t194 - t231 * t197;
t204 = qJD(5) + t206;
t186 = t219 * pkin(4) - t238;
t233 = cos(qJ(6));
t229 = sin(qJ(6));
t203 = qJD(6) + t204;
t200 = t207 * t251 - t230 * t219;
t199 = t230 * t207 + t219 * t251;
t192 = -t229 * t199 + t233 * t200;
t191 = t233 * t199 + t229 * t200;
t184 = t199 * pkin(5) + t186;
t183 = -t199 * pkin(11) + t248;
t182 = t204 * pkin(5) - t200 * pkin(11) + t239;
t1 = [t236 / 0.2e1, 0, 0, t232 ^ 2 * t249 / 0.2e1, t232 * t242, t224 * t241, t224 * t240, t224 ^ 2 / 0.2e1, pkin(1) * t242 + t224 * t237, -pkin(1) * t232 * t249 - t224 * t246, -t201 * t240 + t209 * t215, t202 * t240 + t209 * t216, -t201 * t216 - t202 * t215, t202 ^ 2 / 0.2e1 + t201 ^ 2 / 0.2e1 + t209 ^ 2 / 0.2e1, t207 ^ 2 / 0.2e1, -t207 * t206, -t207 * t219, t206 * t219, t219 ^ 2 / 0.2e1, t205 * t206 - t219 * t238, t205 * t207 + t219 * t247, t200 ^ 2 / 0.2e1, -t200 * t199, t200 * t204, -t199 * t204, t204 ^ 2 / 0.2e1, t186 * t199 + t204 * t239, t186 * t200 - t204 * t248, t192 ^ 2 / 0.2e1, -t192 * t191, t192 * t203, -t191 * t203, t203 ^ 2 / 0.2e1 (t233 * t182 - t229 * t183) * t203 + t184 * t191 -(t229 * t182 + t233 * t183) * t203 + t184 * t192;];
T_reg  = t1;
