% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% T_reg [1x33]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:34:08
% EndTime: 2019-03-09 22:34:08
% DurationCPUTime: 0.18s
% Computational Cost: add. (899->56), mult. (2109->116), div. (0->0), fcn. (1701->12), ass. (0->52)
t256 = cos(qJ(4));
t232 = sin(pkin(6));
t242 = qJD(1) ^ 2;
t255 = t232 ^ 2 * t242;
t250 = cos(pkin(6)) * qJD(1);
t229 = qJD(2) + t250;
t237 = sin(qJ(3));
t240 = cos(qJ(3));
t238 = sin(qJ(2));
t251 = qJD(1) * t232;
t247 = t238 * t251;
t218 = -t240 * t229 + t237 * t247;
t219 = t237 * t229 + t240 * t247;
t236 = sin(qJ(4));
t209 = -t236 * t218 + t256 * t219;
t241 = cos(qJ(2));
t246 = t241 * t251;
t224 = -qJD(3) + t246;
t221 = -qJD(4) + t224;
t249 = pkin(1) * t250;
t252 = pkin(8) * t246 + t238 * t249;
t215 = t229 * pkin(9) + t252;
t217 = (-pkin(2) * t241 - pkin(9) * t238 - pkin(1)) * t251;
t244 = -t237 * t215 + t240 * t217;
t204 = -t224 * pkin(3) - t219 * pkin(10) + t244;
t253 = t240 * t215 + t237 * t217;
t206 = -t218 * pkin(10) + t253;
t245 = t256 * t204 - t236 * t206;
t193 = -t221 * pkin(4) - t209 * qJ(5) + t245;
t208 = t256 * t218 + t236 * t219;
t254 = t236 * t204 + t256 * t206;
t195 = -t208 * qJ(5) + t254;
t231 = sin(pkin(12));
t233 = cos(pkin(12));
t190 = t231 * t193 + t233 * t195;
t248 = t241 * t255;
t199 = -t233 * t208 - t231 * t209;
t189 = t233 * t193 - t231 * t195;
t243 = -pkin(8) * t247 + t241 * t249;
t214 = -t229 * pkin(2) - t243;
t210 = t218 * pkin(3) + t214;
t201 = t208 * pkin(4) + qJD(5) + t210;
t239 = cos(qJ(6));
t235 = sin(qJ(6));
t200 = -t231 * t208 + t233 * t209;
t198 = qJD(6) - t199;
t197 = t239 * t200 - t235 * t221;
t196 = t235 * t200 + t239 * t221;
t191 = -t199 * pkin(5) - t200 * pkin(11) + t201;
t188 = -t221 * pkin(11) + t190;
t187 = t221 * pkin(5) - t189;
t1 = [t242 / 0.2e1, 0, 0, t238 ^ 2 * t255 / 0.2e1, t238 * t248, t229 * t247, t229 * t246, t229 ^ 2 / 0.2e1, pkin(1) * t248 + t243 * t229, -pkin(1) * t238 * t255 - t252 * t229, t219 ^ 2 / 0.2e1, -t219 * t218, -t219 * t224, t218 * t224, t224 ^ 2 / 0.2e1, t214 * t218 - t244 * t224, t214 * t219 + t253 * t224, t209 ^ 2 / 0.2e1, -t209 * t208, -t209 * t221, t208 * t221, t221 ^ 2 / 0.2e1, t210 * t208 - t245 * t221, t210 * t209 + t254 * t221, -t189 * t200 + t190 * t199, t190 ^ 2 / 0.2e1 + t189 ^ 2 / 0.2e1 + t201 ^ 2 / 0.2e1, t197 ^ 2 / 0.2e1, -t197 * t196, t197 * t198, -t196 * t198, t198 ^ 2 / 0.2e1 (-t235 * t188 + t239 * t191) * t198 + t187 * t196 -(t239 * t188 + t235 * t191) * t198 + t187 * t197;];
T_reg  = t1;
