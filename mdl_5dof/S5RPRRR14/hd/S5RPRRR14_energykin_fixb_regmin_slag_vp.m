% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR14_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:19:42
% EndTime: 2019-12-31 19:19:42
% DurationCPUTime: 0.20s
% Computational Cost: add. (489->57), mult. (1627->119), div. (0->0), fcn. (1347->12), ass. (0->50)
t229 = sin(pkin(11));
t254 = cos(pkin(5));
t243 = qJD(1) * t254;
t242 = pkin(1) * t243;
t232 = cos(pkin(11));
t231 = sin(pkin(5));
t248 = qJD(1) * t231;
t245 = t232 * t248;
t221 = qJ(2) * t245 + t229 * t242;
t233 = cos(pkin(6));
t230 = sin(pkin(6));
t244 = t254 * t230;
t251 = t231 * t232;
t209 = (t233 * t251 + t244) * qJD(1) * pkin(8) + t221;
t226 = t232 * t242;
t252 = t229 * t231;
t211 = t226 + (t254 * pkin(2) + (-pkin(8) * t233 - qJ(2)) * t252) * qJD(1);
t236 = sin(qJ(3));
t239 = cos(qJ(3));
t216 = qJD(2) + (-pkin(8) * t229 * t230 - pkin(2) * t232 - pkin(1)) * t248;
t253 = t216 * t230;
t255 = -t236 * t209 + (t211 * t233 + t253) * t239;
t250 = t233 * t236;
t201 = -t230 * t211 + t233 * t216;
t246 = t229 * t248;
t212 = t236 * t246 + (-t230 * t243 - t233 * t245) * t239;
t213 = (t236 * t244 + (t229 * t239 + t232 * t250) * t231) * qJD(1);
t195 = t212 * pkin(3) - t213 * pkin(9) + t201;
t218 = t230 * t245 - t233 * t243 - qJD(3);
t247 = t239 * t209 + t211 * t250 + t236 * t253;
t198 = -t218 * pkin(9) + t247;
t235 = sin(qJ(4));
t238 = cos(qJ(4));
t249 = t235 * t195 + t238 * t198;
t203 = t235 * t213 + t238 * t218;
t241 = t238 * t195 - t235 * t198;
t197 = t218 * pkin(3) - t255;
t237 = cos(qJ(5));
t234 = sin(qJ(5));
t227 = -pkin(1) * t248 + qJD(2);
t220 = -qJ(2) * t246 + t226;
t210 = qJD(4) + t212;
t204 = t238 * t213 - t235 * t218;
t202 = qJD(5) + t203;
t200 = t237 * t204 + t234 * t210;
t199 = t234 * t204 - t237 * t210;
t193 = t203 * pkin(4) - t204 * pkin(10) + t197;
t192 = t210 * pkin(10) + t249;
t191 = -t210 * pkin(4) - t241;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, (t254 * t220 - t227 * t251) * qJD(1), (-t254 * t221 + t227 * t252) * qJD(1), (-t220 * t229 + t221 * t232) * t248, t221 ^ 2 / 0.2e1 + t220 ^ 2 / 0.2e1 + t227 ^ 2 / 0.2e1, t213 ^ 2 / 0.2e1, -t213 * t212, -t213 * t218, t212 * t218, t218 ^ 2 / 0.2e1, t201 * t212 - t255 * t218, t201 * t213 + t247 * t218, t204 ^ 2 / 0.2e1, -t204 * t203, t204 * t210, -t203 * t210, t210 ^ 2 / 0.2e1, t197 * t203 + t241 * t210, t197 * t204 - t249 * t210, t200 ^ 2 / 0.2e1, -t200 * t199, t200 * t202, -t199 * t202, t202 ^ 2 / 0.2e1, (-t234 * t192 + t237 * t193) * t202 + t191 * t199, -(t237 * t192 + t234 * t193) * t202 + t191 * t200;];
T_reg = t1;
