% Calculate minimal parameter regressor of potential energy for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:20:58
% EndTime: 2019-03-09 09:20:58
% DurationCPUTime: 0.15s
% Computational Cost: add. (106->56), mult. (247->86), div. (0->0), fcn. (293->10), ass. (0->35)
t226 = sin(qJ(1));
t245 = g(1) * t226;
t230 = cos(qJ(1));
t244 = g(2) * t230;
t221 = sin(pkin(6));
t225 = sin(qJ(2));
t243 = t221 * t225;
t242 = t221 * t226;
t229 = cos(qJ(2));
t241 = t221 * t229;
t240 = t221 * t230;
t239 = t226 * t225;
t238 = t226 * t229;
t237 = t230 * t225;
t236 = t230 * t229;
t222 = cos(pkin(6));
t235 = pkin(2) * t243 + t222 * pkin(8) + pkin(7);
t234 = -t244 + t245;
t211 = -t222 * t236 + t239;
t212 = t222 * t237 + t238;
t233 = t226 * pkin(1) + t212 * pkin(2) + t211 * qJ(3);
t213 = t222 * t238 + t237;
t214 = -t222 * t239 + t236;
t232 = t230 * pkin(1) + t214 * pkin(2) + pkin(8) * t242 + t213 * qJ(3);
t203 = -g(1) * t213 - g(2) * t211 + g(3) * t241;
t231 = g(1) * t214 + g(2) * t212 + g(3) * t243;
t228 = cos(qJ(5));
t227 = cos(qJ(6));
t224 = sin(qJ(5));
t223 = sin(qJ(6));
t210 = -t222 * t224 - t228 * t241;
t207 = g(3) * t222 + t234 * t221;
t206 = t213 * t228 - t224 * t242;
t205 = t211 * t228 + t224 * t240;
t1 = [0, -g(1) * t230 - g(2) * t226, t234, 0, 0, 0, 0, 0, -t231, -t203, -t231, -t207, t203, -g(1) * t232 - g(2) * (-pkin(8) * t240 + t233) - g(3) * (-qJ(3) * t241 + t235) t203, t231, t207, -g(1) * (t214 * pkin(3) + t232) - g(2) * (t212 * pkin(3) + t233) - g(3) * (-t222 * qJ(4) + t235) + (qJ(4) * t245 - g(3) * (pkin(3) * t225 - qJ(3) * t229) - (-pkin(8) + qJ(4)) * t244) * t221, 0, 0, 0, 0, 0, -g(1) * t206 - g(2) * t205 - g(3) * t210, -g(1) * (-t213 * t224 - t228 * t242) - g(2) * (-t211 * t224 + t228 * t240) - g(3) * (-t222 * t228 + t224 * t241) 0, 0, 0, 0, 0, -g(1) * (t206 * t227 + t214 * t223) - g(2) * (t205 * t227 + t212 * t223) - g(3) * (t210 * t227 + t223 * t243) -g(1) * (-t206 * t223 + t214 * t227) - g(2) * (-t205 * t223 + t212 * t227) - g(3) * (-t210 * t223 + t227 * t243);];
U_reg  = t1;
