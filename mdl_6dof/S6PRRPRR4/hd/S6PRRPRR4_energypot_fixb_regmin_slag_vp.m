% Calculate minimal parameter regressor of potential energy for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:14:27
% EndTime: 2019-03-08 22:14:27
% DurationCPUTime: 0.13s
% Computational Cost: add. (142->58), mult. (351->97), div. (0->0), fcn. (449->12), ass. (0->36)
t226 = sin(pkin(6));
t245 = pkin(7) * t226;
t231 = sin(qJ(3));
t244 = t226 * t231;
t232 = sin(qJ(2));
t243 = t226 * t232;
t235 = cos(qJ(3));
t242 = t226 * t235;
t236 = cos(qJ(2));
t241 = t226 * t236;
t228 = cos(pkin(6));
t240 = t228 * t232;
t239 = t228 * t236;
t225 = sin(pkin(11));
t227 = cos(pkin(11));
t217 = t225 * t236 + t227 * t240;
t212 = t217 * t231 + t227 * t242;
t219 = -t225 * t240 + t227 * t236;
t214 = t219 * t231 - t225 * t242;
t220 = -t228 * t235 + t231 * t243;
t238 = g(1) * t214 + g(2) * t212 + g(3) * t220;
t216 = t225 * t232 - t227 * t239;
t218 = t225 * t239 + t227 * t232;
t237 = -g(1) * t218 - g(2) * t216 + g(3) * t241;
t234 = cos(qJ(5));
t233 = cos(qJ(6));
t230 = sin(qJ(5));
t229 = sin(qJ(6));
t221 = t228 * t231 + t232 * t242;
t215 = t219 * t235 + t225 * t244;
t213 = t217 * t235 - t227 * t244;
t211 = t220 * t230 + t221 * t234;
t210 = t214 * t230 + t215 * t234;
t209 = t212 * t230 + t213 * t234;
t208 = -g(1) * t215 - g(2) * t213 - g(3) * t221;
t1 = [-g(3) * qJ(1), 0, -g(1) * t219 - g(2) * t217 - g(3) * t243, -t237, 0, 0, 0, 0, 0, t208, t238, t208, t237, -t238, -g(1) * (t227 * pkin(1) + t219 * pkin(2) + t215 * pkin(3) + t218 * pkin(8) + t214 * qJ(4) + t225 * t245) - g(2) * (t225 * pkin(1) + t217 * pkin(2) + t213 * pkin(3) + t216 * pkin(8) + t212 * qJ(4) - t227 * t245) - g(3) * (t221 * pkin(3) + t228 * pkin(7) + t220 * qJ(4) + qJ(1) + (pkin(2) * t232 - pkin(8) * t236) * t226) 0, 0, 0, 0, 0, -g(1) * t210 - g(2) * t209 - g(3) * t211, -g(1) * (t214 * t234 - t215 * t230) - g(2) * (t212 * t234 - t213 * t230) - g(3) * (t220 * t234 - t221 * t230) 0, 0, 0, 0, 0, -g(1) * (t210 * t233 - t218 * t229) - g(2) * (t209 * t233 - t216 * t229) - g(3) * (t211 * t233 + t229 * t241) -g(1) * (-t210 * t229 - t218 * t233) - g(2) * (-t209 * t229 - t216 * t233) - g(3) * (-t211 * t229 + t233 * t241);];
U_reg  = t1;
