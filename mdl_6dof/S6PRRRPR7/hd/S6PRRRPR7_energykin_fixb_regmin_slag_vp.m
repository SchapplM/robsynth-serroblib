% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRPR7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_energykin_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:44:36
% EndTime: 2019-03-08 23:44:37
% DurationCPUTime: 0.22s
% Computational Cost: add. (548->52), mult. (1304->113), div. (0->0), fcn. (1063->14), ass. (0->50)
t244 = cos(qJ(2));
t256 = qJD(1) * sin(pkin(6));
t224 = qJD(2) * pkin(2) + t244 * t256;
t233 = sin(pkin(7));
t235 = cos(pkin(7));
t255 = qJD(1) * cos(pkin(6));
t246 = t224 * t235 + t233 * t255;
t260 = cos(pkin(13));
t245 = qJD(2) ^ 2;
t258 = t233 ^ 2 * t245;
t243 = cos(qJ(3));
t254 = qJD(2) * t233;
t249 = t243 * t254;
t227 = -qJD(4) + t249;
t230 = qJD(2) * t235 + qJD(3);
t240 = sin(qJ(2));
t222 = pkin(9) * t254 + t240 * t256;
t239 = sin(qJ(3));
t253 = t243 * t222 + t239 * t246;
t209 = pkin(10) * t230 + t253;
t229 = t235 * t255;
t212 = t229 + (-t224 + (-pkin(3) * t243 - pkin(10) * t239) * qJD(2)) * t233;
t238 = sin(qJ(4));
t242 = cos(qJ(4));
t257 = t209 * t242 + t212 * t238;
t200 = -qJ(5) * t227 + t257;
t220 = t239 * t222;
t208 = -t230 * pkin(3) - t243 * t246 + t220;
t251 = t239 * t254;
t217 = -t230 * t242 + t238 * t251;
t218 = t230 * t238 + t242 * t251;
t203 = t217 * pkin(4) - t218 * qJ(5) + t208;
t232 = sin(pkin(13));
t196 = t200 * t260 + t203 * t232;
t250 = (-t233 * t224 + t229) * t254;
t248 = qJD(2) * t256;
t195 = -t200 * t232 + t203 * t260;
t247 = -t209 * t238 + t242 * t212;
t199 = pkin(4) * t227 + qJD(5) - t247;
t241 = cos(qJ(6));
t237 = sin(qJ(6));
t216 = qJD(6) + t217;
t214 = t218 * t260 - t227 * t232;
t213 = t218 * t232 + t227 * t260;
t205 = -t213 * t237 + t214 * t241;
t204 = t213 * t241 + t214 * t237;
t197 = pkin(5) * t213 + t199;
t194 = -pkin(11) * t213 + t196;
t193 = pkin(5) * t217 - pkin(11) * t214 + t195;
t1 = [qJD(1) ^ 2 / 0.2e1, t245 / 0.2e1, t244 * t248, -t240 * t248, t239 ^ 2 * t258 / 0.2e1, t239 * t243 * t258, t230 * t251, t230 * t249, t230 ^ 2 / 0.2e1, -t220 * t230 + (t230 * t246 - t250) * t243, -t230 * t253 + t239 * t250, t218 ^ 2 / 0.2e1, -t218 * t217, -t218 * t227, t217 * t227, t227 ^ 2 / 0.2e1, t208 * t217 - t227 * t247, t208 * t218 + t227 * t257, t195 * t217 + t199 * t213, -t196 * t217 + t199 * t214, -t195 * t214 - t196 * t213, t196 ^ 2 / 0.2e1 + t195 ^ 2 / 0.2e1 + t199 ^ 2 / 0.2e1, t205 ^ 2 / 0.2e1, -t205 * t204, t205 * t216, -t204 * t216, t216 ^ 2 / 0.2e1 (t193 * t241 - t194 * t237) * t216 + t197 * t204 -(t193 * t237 + t194 * t241) * t216 + t197 * t205;];
T_reg  = t1;
