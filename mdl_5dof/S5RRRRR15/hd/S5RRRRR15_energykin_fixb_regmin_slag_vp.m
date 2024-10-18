% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% T_reg [1x31]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRRR15_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:18
% EndTime: 2024-09-27 22:26:18
% DurationCPUTime: 0.00s
% Computational Cost: add. (481->50), mult. (1426->107), div. (0->0), fcn. (1177->12), ass. (0->49)
t231 = sin(qJ(3));
t234 = cos(qJ(3));
t235 = cos(qJ(2));
t226 = sin(pkin(5));
t246 = qJD(1) * t226;
t241 = t235 * t246;
t232 = sin(qJ(2));
t242 = t232 * t246;
t212 = t231 * t242 - t234 * t241;
t213 = (t231 * t235 + t232 * t234) * t246;
t230 = sin(qJ(4));
t254 = cos(qJ(4));
t203 = t254 * t212 + t230 * t213;
t245 = cos(pkin(5)) * qJD(1);
t223 = qJD(2) + t245;
t219 = qJD(3) + t223;
t218 = qJD(4) + t219;
t225 = sin(pkin(6));
t227 = cos(pkin(6));
t237 = -t203 * t227 + t218 * t225;
t204 = -t230 * t212 + t254 * t213;
t253 = pkin(11) * t204;
t236 = qJD(1) ^ 2;
t250 = t226 ^ 2 * t236;
t244 = pkin(1) * t245;
t222 = t235 * t244;
t208 = t223 * pkin(2) + t222 + (-pkin(8) - pkin(9)) * t242;
t247 = pkin(8) * t241 + t232 * t244;
t210 = pkin(9) * t241 + t247;
t239 = t234 * t208 - t231 * t210;
t197 = t219 * pkin(3) - t213 * pkin(10) + t239;
t248 = t231 * t208 + t234 * t210;
t199 = -t212 * pkin(10) + t248;
t249 = t230 * t197 + t254 * t199;
t243 = t235 * t250;
t240 = t254 * t197 - t230 * t199;
t191 = t218 * pkin(4) - t227 * t253 + t240;
t215 = (-pkin(2) * t235 - pkin(1)) * t246;
t207 = t212 * pkin(3) + t215;
t192 = t203 * pkin(4) - t225 * t253 + t207;
t238 = t191 * t227 + t192 * t225;
t233 = cos(qJ(5));
t229 = sin(qJ(5));
t200 = -t225 * t203 - t227 * t218 - qJD(5);
t194 = t233 * t204 + t237 * t229;
t193 = t229 * t204 - t237 * t233;
t190 = t237 * pkin(11) + t249;
t189 = -t225 * t191 + t227 * t192;
t1 = [t236 / 0.2e1, 0, 0, t232 ^ 2 * t250 / 0.2e1, t232 * t243, t223 * t242, t223 * t241, t223 ^ 2 / 0.2e1, pkin(1) * t243 + (-pkin(8) * t242 + t222) * t223, -pkin(1) * t232 * t250 - t247 * t223, t213 ^ 2 / 0.2e1, -t213 * t212, t213 * t219, -t212 * t219, t219 ^ 2 / 0.2e1, t215 * t212 + t239 * t219, t215 * t213 - t248 * t219, t204 ^ 2 / 0.2e1, -t204 * t203, t204 * t218, -t203 * t218, t218 ^ 2 / 0.2e1, t207 * t203 + t240 * t218, t207 * t204 - t249 * t218, t194 ^ 2 / 0.2e1, -t194 * t193, -t194 * t200, t193 * t200, t200 ^ 2 / 0.2e1, -(-t229 * t190 + t238 * t233) * t200 + t189 * t193, (t233 * t190 + t238 * t229) * t200 + t189 * t194;];
T_reg = t1;
