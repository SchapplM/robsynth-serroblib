% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% T_reg [1x38]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRRRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:04:46
% EndTime: 2019-03-10 04:04:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (872->57), mult. (2005->119), div. (0->0), fcn. (1655->12), ass. (0->53)
t257 = cos(qJ(3));
t256 = cos(qJ(4));
t230 = sin(pkin(6));
t240 = qJD(1) ^ 2;
t255 = t230 ^ 2 * t240;
t249 = cos(pkin(6)) * qJD(1);
t228 = qJD(2) + t249;
t235 = sin(qJ(3));
t236 = sin(qJ(2));
t250 = qJD(1) * t230;
t246 = t236 * t250;
t215 = -t228 * t257 + t235 * t246;
t216 = t235 * t228 + t246 * t257;
t234 = sin(qJ(4));
t206 = -t234 * t215 + t216 * t256;
t239 = cos(qJ(2));
t245 = t239 * t250;
t223 = -qJD(3) + t245;
t220 = -qJD(4) + t223;
t248 = pkin(1) * t249;
t251 = pkin(8) * t245 + t236 * t248;
t212 = t228 * pkin(9) + t251;
t214 = (-pkin(2) * t239 - pkin(9) * t236 - pkin(1)) * t250;
t243 = -t235 * t212 + t257 * t214;
t201 = -t223 * pkin(3) - t216 * pkin(10) + t243;
t252 = t212 * t257 + t235 * t214;
t203 = -t215 * pkin(10) + t252;
t244 = t256 * t201 - t234 * t203;
t190 = -t220 * pkin(4) - t206 * pkin(11) + t244;
t205 = t215 * t256 + t234 * t216;
t253 = t234 * t201 + t256 * t203;
t192 = -t205 * pkin(11) + t253;
t233 = sin(qJ(5));
t238 = cos(qJ(5));
t254 = t233 * t190 + t238 * t192;
t247 = t239 * t255;
t196 = t238 * t205 + t233 * t206;
t242 = t238 * t190 - t233 * t192;
t241 = -pkin(8) * t246 + t239 * t248;
t211 = -t228 * pkin(2) - t241;
t207 = t215 * pkin(3) + t211;
t198 = t205 * pkin(4) + t207;
t237 = cos(qJ(6));
t232 = sin(qJ(6));
t218 = -qJD(5) + t220;
t197 = -t233 * t205 + t238 * t206;
t195 = qJD(6) + t196;
t194 = t237 * t197 - t232 * t218;
t193 = t232 * t197 + t237 * t218;
t188 = t196 * pkin(5) - t197 * pkin(12) + t198;
t187 = -t218 * pkin(12) + t254;
t186 = t218 * pkin(5) - t242;
t1 = [t240 / 0.2e1, 0, 0, t236 ^ 2 * t255 / 0.2e1, t236 * t247, t228 * t246, t228 * t245, t228 ^ 2 / 0.2e1, pkin(1) * t247 + t228 * t241, -pkin(1) * t236 * t255 - t228 * t251, t216 ^ 2 / 0.2e1, -t216 * t215, -t216 * t223, t215 * t223, t223 ^ 2 / 0.2e1, t211 * t215 - t223 * t243, t211 * t216 + t223 * t252, t206 ^ 2 / 0.2e1, -t206 * t205, -t206 * t220, t205 * t220, t220 ^ 2 / 0.2e1, t207 * t205 - t220 * t244, t207 * t206 + t220 * t253, t197 ^ 2 / 0.2e1, -t197 * t196, -t197 * t218, t196 * t218, t218 ^ 2 / 0.2e1, t198 * t196 - t218 * t242, t198 * t197 + t218 * t254, t194 ^ 2 / 0.2e1, -t194 * t193, t194 * t195, -t193 * t195, t195 ^ 2 / 0.2e1 (-t232 * t187 + t237 * t188) * t195 + t186 * t193 -(t237 * t187 + t232 * t188) * t195 + t186 * t194;];
T_reg  = t1;
