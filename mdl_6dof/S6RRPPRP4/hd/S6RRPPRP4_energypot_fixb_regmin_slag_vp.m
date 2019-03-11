% Calculate minimal parameter regressor of potential energy for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:39:41
% EndTime: 2019-03-09 08:39:41
% DurationCPUTime: 0.12s
% Computational Cost: add. (142->61), mult. (313->83), div. (0->0), fcn. (364->8), ass. (0->39)
t226 = cos(pkin(9));
t228 = sin(qJ(2));
t246 = t226 * t228;
t225 = sin(pkin(9));
t247 = t225 * t228;
t249 = pkin(3) * t246 + qJ(4) * t247;
t248 = t228 * pkin(2) + pkin(6);
t229 = sin(qJ(1));
t245 = t228 * t229;
t232 = cos(qJ(1));
t244 = t228 * t232;
t231 = cos(qJ(2));
t243 = t229 * t231;
t242 = t232 * t225;
t241 = t232 * t226;
t240 = t229 * pkin(7) + qJ(3) * t244 + (pkin(2) * t231 + pkin(1)) * t232;
t239 = -t231 * qJ(3) + t248;
t238 = g(1) * t232 + g(2) * t229;
t237 = t229 * pkin(1) + pkin(2) * t243 - t232 * pkin(7) + qJ(3) * t245;
t208 = t225 * t243 + t241;
t209 = t226 * t243 - t242;
t227 = sin(qJ(5));
t230 = cos(qJ(5));
t197 = -t208 * t230 + t209 * t227;
t210 = -t229 * t226 + t231 * t242;
t211 = t229 * t225 + t231 * t241;
t199 = -t210 * t230 + t211 * t227;
t204 = t227 * t246 - t230 * t247;
t236 = g(1) * t199 + g(2) * t197 + g(3) * t204;
t235 = t211 * pkin(3) + t210 * qJ(4) + t240;
t234 = g(1) * t210 + g(2) * t208 + g(3) * t247;
t233 = t209 * pkin(3) + t208 * qJ(4) + t237;
t205 = (t225 * t227 + t226 * t230) * t228;
t203 = -g(3) * t231 + t238 * t228;
t200 = t210 * t227 + t211 * t230;
t198 = t208 * t227 + t209 * t230;
t196 = -g(1) * t211 - g(2) * t209 - g(3) * t246;
t195 = -g(1) * t200 - g(2) * t198 - g(3) * t205;
t1 = [0, -t238, g(1) * t229 - g(2) * t232, 0, 0, 0, 0, 0, -g(3) * t228 - t238 * t231, t203, t196, t234, -t203, -g(1) * t240 - g(2) * t237 - g(3) * t239, t196, -t203, -t234, -g(1) * t235 - g(2) * t233 - g(3) * (t239 + t249) 0, 0, 0, 0, 0, t195, t236, t195, t203, -t236, -g(1) * (t211 * pkin(4) + t200 * pkin(5) - pkin(8) * t244 + t199 * qJ(6) + t235) - g(2) * (t209 * pkin(4) + t198 * pkin(5) - pkin(8) * t245 + t197 * qJ(6) + t233) - g(3) * (pkin(4) * t246 + t205 * pkin(5) + t204 * qJ(6) + (pkin(8) - qJ(3)) * t231 + t248 + t249);];
U_reg  = t1;
