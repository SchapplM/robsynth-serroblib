% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% T_reg [1x35]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRPRR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:46:28
% EndTime: 2019-03-09 19:46:28
% DurationCPUTime: 0.17s
% Computational Cost: add. (869->58), mult. (1993->120), div. (0->0), fcn. (1611->12), ass. (0->52)
t249 = cos(qJ(5));
t225 = sin(pkin(6));
t235 = qJD(1) ^ 2;
t248 = t225 ^ 2 * t235;
t243 = cos(pkin(6)) * qJD(1);
t222 = qJD(2) + t243;
t234 = cos(qJ(2));
t231 = sin(qJ(2));
t244 = qJD(1) * t225;
t240 = t231 * t244;
t242 = pkin(1) * t243;
t236 = -pkin(8) * t240 + t234 * t242;
t208 = -t222 * pkin(2) - t236;
t230 = sin(qJ(3));
t233 = cos(qJ(3));
t213 = -t233 * t222 + t230 * t240;
t214 = t230 * t222 + t233 * t240;
t198 = t213 * pkin(3) - t214 * qJ(4) + t208;
t239 = t234 * t244;
t217 = -qJD(3) + t239;
t245 = pkin(8) * t239 + t231 * t242;
t209 = t222 * pkin(9) + t245;
t211 = (-pkin(2) * t234 - pkin(9) * t231 - pkin(1)) * t244;
t246 = t233 * t209 + t230 * t211;
t201 = -t217 * qJ(4) + t246;
t224 = sin(pkin(12));
t226 = cos(pkin(12));
t190 = t226 * t198 - t224 * t201;
t204 = t226 * t214 - t224 * t217;
t184 = t213 * pkin(4) - t204 * pkin(10) + t190;
t191 = t224 * t198 + t226 * t201;
t203 = t224 * t214 + t226 * t217;
t187 = -t203 * pkin(10) + t191;
t229 = sin(qJ(5));
t247 = t229 * t184 + t249 * t187;
t241 = t234 * t248;
t238 = t249 * t184 - t229 * t187;
t237 = -t230 * t209 + t233 * t211;
t212 = qJD(5) + t213;
t200 = t217 * pkin(3) + qJD(4) - t237;
t192 = t203 * pkin(4) + t200;
t232 = cos(qJ(6));
t228 = sin(qJ(6));
t210 = qJD(6) + t212;
t195 = -t229 * t203 + t249 * t204;
t194 = t249 * t203 + t229 * t204;
t189 = -t228 * t194 + t232 * t195;
t188 = t232 * t194 + t228 * t195;
t185 = t194 * pkin(5) + t192;
t181 = -t194 * pkin(11) + t247;
t180 = t212 * pkin(5) - t195 * pkin(11) + t238;
t1 = [t235 / 0.2e1, 0, 0, t231 ^ 2 * t248 / 0.2e1, t231 * t241, t222 * t240, t222 * t239, t222 ^ 2 / 0.2e1, pkin(1) * t241 + t236 * t222, -pkin(1) * t231 * t248 - t245 * t222, t214 ^ 2 / 0.2e1, -t214 * t213, -t214 * t217, t213 * t217, t217 ^ 2 / 0.2e1, t208 * t213 - t237 * t217, t208 * t214 + t246 * t217, t190 * t213 + t200 * t203, -t191 * t213 + t200 * t204, -t190 * t204 - t191 * t203, t191 ^ 2 / 0.2e1 + t190 ^ 2 / 0.2e1 + t200 ^ 2 / 0.2e1, t195 ^ 2 / 0.2e1, -t195 * t194, t195 * t212, -t194 * t212, t212 ^ 2 / 0.2e1, t192 * t194 + t212 * t238, t192 * t195 - t247 * t212, t189 ^ 2 / 0.2e1, -t189 * t188, t189 * t210, -t188 * t210, t210 ^ 2 / 0.2e1 (t232 * t180 - t228 * t181) * t210 + t185 * t188 -(t228 * t180 + t232 * t181) * t210 + t185 * t189;];
T_reg  = t1;
