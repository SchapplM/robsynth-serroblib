% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% T_reg [1x29]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6PRRRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:34:24
% EndTime: 2019-03-09 00:34:24
% DurationCPUTime: 0.17s
% Computational Cost: add. (477->48), mult. (1125->105), div. (0->0), fcn. (904->12), ass. (0->46)
t234 = cos(qJ(2));
t247 = qJD(1) * sin(pkin(6));
t215 = qJD(2) * pkin(2) + t234 * t247;
t223 = sin(pkin(7));
t225 = cos(pkin(7));
t246 = qJD(1) * cos(pkin(6));
t236 = t215 * t225 + t223 * t246;
t235 = qJD(2) ^ 2;
t250 = t223 ^ 2 * t235;
t233 = cos(qJ(3));
t245 = qJD(2) * t223;
t240 = t233 * t245;
t218 = -qJD(4) + t240;
t221 = t225 * qJD(2) + qJD(3);
t230 = sin(qJ(2));
t213 = pkin(9) * t245 + t230 * t247;
t229 = sin(qJ(3));
t244 = t233 * t213 + t236 * t229;
t201 = t221 * pkin(10) + t244;
t220 = t225 * t246;
t203 = t220 + (-t215 + (-pkin(3) * t233 - pkin(10) * t229) * qJD(2)) * t223;
t228 = sin(qJ(4));
t232 = cos(qJ(4));
t248 = t232 * t201 + t228 * t203;
t195 = -t218 * pkin(11) + t248;
t211 = t229 * t213;
t200 = -t221 * pkin(3) - t236 * t233 + t211;
t242 = t229 * t245;
t208 = -t232 * t221 + t228 * t242;
t209 = t228 * t221 + t232 * t242;
t197 = t208 * pkin(4) - t209 * pkin(11) + t200;
t227 = sin(qJ(5));
t231 = cos(qJ(5));
t249 = t231 * t195 + t227 * t197;
t241 = (-t223 * t215 + t220) * t245;
t239 = qJD(2) * t247;
t238 = -t228 * t201 + t232 * t203;
t237 = -t227 * t195 + t231 * t197;
t194 = t218 * pkin(4) - t238;
t207 = qJD(5) + t208;
t205 = t231 * t209 - t227 * t218;
t204 = t227 * t209 + t231 * t218;
t192 = t204 * pkin(5) - t205 * qJ(6) + t194;
t191 = t207 * qJ(6) + t249;
t190 = -t207 * pkin(5) + qJD(6) - t237;
t1 = [qJD(1) ^ 2 / 0.2e1, t235 / 0.2e1, t234 * t239, -t230 * t239, t229 ^ 2 * t250 / 0.2e1, t229 * t233 * t250, t221 * t242, t221 * t240, t221 ^ 2 / 0.2e1, -t211 * t221 + (t236 * t221 - t241) * t233, -t244 * t221 + t229 * t241, t209 ^ 2 / 0.2e1, -t209 * t208, -t209 * t218, t208 * t218, t218 ^ 2 / 0.2e1, t200 * t208 - t238 * t218, t200 * t209 + t248 * t218, t205 ^ 2 / 0.2e1, -t205 * t204, t205 * t207, -t204 * t207, t207 ^ 2 / 0.2e1, t194 * t204 + t237 * t207, t194 * t205 - t249 * t207, -t190 * t207 + t192 * t204, t190 * t205 - t191 * t204, t191 * t207 - t192 * t205, t191 ^ 2 / 0.2e1 + t192 ^ 2 / 0.2e1 + t190 ^ 2 / 0.2e1;];
T_reg  = t1;
