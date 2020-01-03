% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:55:11
% EndTime: 2019-12-31 22:55:11
% DurationCPUTime: 0.22s
% Computational Cost: add. (592->52), mult. (1611->112), div. (0->0), fcn. (1336->12), ass. (0->49)
t248 = cos(pkin(5)) * qJD(1);
t225 = qJD(2) + t248;
t227 = sin(pkin(6));
t229 = cos(pkin(6));
t238 = cos(qJ(2));
t228 = sin(pkin(5));
t249 = qJD(1) * t228;
t243 = t238 * t249;
t257 = t225 * t227 + t229 * t243;
t234 = sin(qJ(2));
t247 = pkin(1) * t248;
t250 = pkin(8) * t243 + t234 * t247;
t208 = t257 * pkin(9) + t250;
t224 = t238 * t247;
t244 = t234 * t249;
t210 = t225 * pkin(2) + t224 + (-pkin(9) * t229 - pkin(8)) * t244;
t216 = (-pkin(9) * t227 * t234 - pkin(2) * t238 - pkin(1)) * t249;
t233 = sin(qJ(3));
t237 = cos(qJ(3));
t256 = -t233 * t208 + (t210 * t229 + t216 * t227) * t237;
t239 = qJD(1) ^ 2;
t254 = t228 ^ 2 * t239;
t253 = t227 * t233;
t252 = t229 * t233;
t200 = -t227 * t210 + t229 * t216;
t211 = t233 * t244 - t257 * t237;
t212 = t225 * t253 + (t234 * t237 + t238 * t252) * t249;
t194 = t211 * pkin(3) - t212 * pkin(10) + t200;
t217 = -t229 * t225 + t227 * t243 - qJD(3);
t245 = t237 * t208 + t210 * t252 + t216 * t253;
t197 = -t217 * pkin(10) + t245;
t232 = sin(qJ(4));
t236 = cos(qJ(4));
t251 = t232 * t194 + t236 * t197;
t246 = t238 * t254;
t202 = t232 * t212 + t236 * t217;
t241 = t236 * t194 - t232 * t197;
t196 = t217 * pkin(3) - t256;
t235 = cos(qJ(5));
t231 = sin(qJ(5));
t209 = qJD(4) + t211;
t203 = t236 * t212 - t232 * t217;
t201 = qJD(5) + t202;
t199 = t235 * t203 + t231 * t209;
t198 = t231 * t203 - t235 * t209;
t192 = t202 * pkin(4) - t203 * pkin(11) + t196;
t191 = t209 * pkin(11) + t251;
t190 = -t209 * pkin(4) - t241;
t1 = [t239 / 0.2e1, 0, 0, t234 ^ 2 * t254 / 0.2e1, t234 * t246, t225 * t244, t225 * t243, t225 ^ 2 / 0.2e1, pkin(1) * t246 + (-pkin(8) * t244 + t224) * t225, -pkin(1) * t234 * t254 - t250 * t225, t212 ^ 2 / 0.2e1, -t212 * t211, -t212 * t217, t211 * t217, t217 ^ 2 / 0.2e1, t200 * t211 - t256 * t217, t200 * t212 + t245 * t217, t203 ^ 2 / 0.2e1, -t203 * t202, t203 * t209, -t202 * t209, t209 ^ 2 / 0.2e1, t196 * t202 + t241 * t209, t196 * t203 - t251 * t209, t199 ^ 2 / 0.2e1, -t199 * t198, t199 * t201, -t198 * t201, t201 ^ 2 / 0.2e1, (-t231 * t191 + t235 * t192) * t201 + t190 * t198, -(t235 * t191 + t231 * t192) * t201 + t190 * t199;];
T_reg = t1;
