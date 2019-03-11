% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:42:54
% EndTime: 2019-03-09 10:42:54
% DurationCPUTime: 0.13s
% Computational Cost: add. (524->55), mult. (1458->110), div. (0->0), fcn. (1134->10), ass. (0->48)
t241 = pkin(4) + pkin(10);
t219 = sin(pkin(6));
t228 = qJD(1) ^ 2;
t240 = t219 ^ 2 * t228;
t227 = cos(qJ(2));
t236 = cos(pkin(6)) * qJD(1);
t235 = pkin(1) * t236;
t215 = t227 * t235;
t216 = qJD(2) + t236;
t224 = sin(qJ(2));
t237 = qJD(1) * t219;
t233 = t224 * t237;
t203 = t216 * pkin(2) + t215 + (-pkin(8) - qJ(3)) * t233;
t232 = t227 * t237;
t238 = pkin(8) * t232 + t224 * t235;
t206 = qJ(3) * t232 + t238;
t218 = sin(pkin(11));
t220 = cos(pkin(11));
t193 = t218 * t203 + t220 * t206;
t191 = t216 * pkin(9) + t193;
t208 = -t218 * t233 + t220 * t232;
t209 = (t218 * t227 + t220 * t224) * t237;
t210 = qJD(3) + (-pkin(2) * t227 - pkin(1)) * t237;
t197 = -t208 * pkin(3) - t209 * pkin(9) + t210;
t223 = sin(qJ(4));
t226 = cos(qJ(4));
t239 = t226 * t191 + t223 * t197;
t234 = t227 * t240;
t231 = -t223 * t191 + t226 * t197;
t192 = t220 * t203 - t218 * t206;
t207 = qJD(4) - t208;
t186 = -t207 * qJ(5) - t239;
t230 = qJD(5) - t231;
t202 = t226 * t209 + t223 * t216;
t190 = -t216 * pkin(3) - t192;
t229 = -t202 * qJ(5) + t190;
t225 = cos(qJ(6));
t222 = sin(qJ(6));
t201 = t223 * t209 - t226 * t216;
t200 = qJD(6) + t202;
t195 = t222 * t201 + t225 * t207;
t194 = -t225 * t201 + t222 * t207;
t187 = t201 * pkin(4) + t229;
t185 = -t207 * pkin(4) + t230;
t184 = t241 * t201 + t229;
t183 = -t201 * pkin(5) - t186;
t182 = t202 * pkin(5) - t241 * t207 + t230;
t1 = [t228 / 0.2e1, 0, 0, t224 ^ 2 * t240 / 0.2e1, t224 * t234, t216 * t233, t216 * t232, t216 ^ 2 / 0.2e1, pkin(1) * t234 + (-pkin(8) * t233 + t215) * t216, -pkin(1) * t224 * t240 - t238 * t216, -t192 * t209 + t193 * t208, t193 ^ 2 / 0.2e1 + t192 ^ 2 / 0.2e1 + t210 ^ 2 / 0.2e1, t202 ^ 2 / 0.2e1, -t202 * t201, t202 * t207, -t201 * t207, t207 ^ 2 / 0.2e1, t190 * t201 + t231 * t207, t190 * t202 - t239 * t207, t185 * t202 + t186 * t201, t185 * t207 - t187 * t201, -t186 * t207 - t187 * t202, t187 ^ 2 / 0.2e1 + t186 ^ 2 / 0.2e1 + t185 ^ 2 / 0.2e1, t195 ^ 2 / 0.2e1, -t195 * t194, t195 * t200, -t194 * t200, t200 ^ 2 / 0.2e1 (t225 * t182 - t222 * t184) * t200 + t183 * t194 -(t222 * t182 + t225 * t184) * t200 + t183 * t195;];
T_reg  = t1;
