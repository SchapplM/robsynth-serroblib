% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RRPRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:05:24
% EndTime: 2019-03-09 12:05:25
% DurationCPUTime: 0.14s
% Computational Cost: add. (531->52), mult. (1476->106), div. (0->0), fcn. (1160->10), ass. (0->47)
t239 = cos(qJ(5));
t218 = sin(pkin(6));
t226 = qJD(1) ^ 2;
t238 = t218 ^ 2 * t226;
t217 = sin(pkin(11));
t219 = cos(pkin(11));
t225 = cos(qJ(2));
t234 = qJD(1) * t218;
t229 = t225 * t234;
t223 = sin(qJ(2));
t230 = t223 * t234;
t207 = -t217 * t230 + t219 * t229;
t206 = qJD(4) - t207;
t233 = cos(pkin(6)) * qJD(1);
t232 = pkin(1) * t233;
t214 = t225 * t232;
t215 = qJD(2) + t233;
t201 = t215 * pkin(2) + t214 + (-pkin(8) - qJ(3)) * t230;
t235 = pkin(8) * t229 + t223 * t232;
t204 = qJ(3) * t229 + t235;
t192 = t217 * t201 + t219 * t204;
t190 = t215 * pkin(9) + t192;
t208 = (t217 * t225 + t219 * t223) * t234;
t209 = qJD(3) + (-pkin(2) * t225 - pkin(1)) * t234;
t196 = -t207 * pkin(3) - t208 * pkin(9) + t209;
t222 = sin(qJ(4));
t224 = cos(qJ(4));
t236 = t224 * t190 + t222 * t196;
t183 = t206 * pkin(10) + t236;
t191 = t219 * t201 - t217 * t204;
t189 = -t215 * pkin(3) - t191;
t199 = t222 * t208 - t224 * t215;
t200 = t224 * t208 + t222 * t215;
t186 = t199 * pkin(4) - t200 * pkin(10) + t189;
t221 = sin(qJ(5));
t237 = t239 * t183 + t221 * t186;
t231 = t225 * t238;
t228 = -t221 * t183 + t239 * t186;
t227 = -t222 * t190 + t224 * t196;
t182 = -t206 * pkin(4) - t227;
t198 = qJD(5) + t199;
t194 = t239 * t200 + t221 * t206;
t193 = t221 * t200 - t239 * t206;
t180 = t193 * pkin(5) + qJD(6) + t182;
t179 = -t193 * qJ(6) + t237;
t178 = t198 * pkin(5) - t194 * qJ(6) + t228;
t1 = [t226 / 0.2e1, 0, 0, t223 ^ 2 * t238 / 0.2e1, t223 * t231, t215 * t230, t215 * t229, t215 ^ 2 / 0.2e1, pkin(1) * t231 + (-pkin(8) * t230 + t214) * t215, -pkin(1) * t223 * t238 - t235 * t215, -t191 * t208 + t192 * t207, t192 ^ 2 / 0.2e1 + t191 ^ 2 / 0.2e1 + t209 ^ 2 / 0.2e1, t200 ^ 2 / 0.2e1, -t200 * t199, t200 * t206, -t199 * t206, t206 ^ 2 / 0.2e1, t189 * t199 + t227 * t206, t189 * t200 - t236 * t206, t194 ^ 2 / 0.2e1, -t194 * t193, t194 * t198, -t193 * t198, t198 ^ 2 / 0.2e1, t182 * t193 + t228 * t198, t182 * t194 - t237 * t198, -t178 * t194 - t179 * t193, t179 ^ 2 / 0.2e1 + t178 ^ 2 / 0.2e1 + t180 ^ 2 / 0.2e1;];
T_reg  = t1;
