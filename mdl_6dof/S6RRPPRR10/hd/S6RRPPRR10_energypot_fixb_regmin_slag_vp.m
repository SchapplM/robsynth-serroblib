% Calculate minimal parameter regressor of potential energy for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:37:00
% EndTime: 2019-03-09 09:37:00
% DurationCPUTime: 0.17s
% Computational Cost: add. (92->52), mult. (132->72), div. (0->0), fcn. (138->10), ass. (0->36)
t227 = sin(qJ(2));
t251 = qJ(3) * t227 + pkin(1);
t229 = cos(qJ(2));
t250 = g(3) * t229;
t224 = pkin(10) + qJ(5);
t219 = qJ(6) + t224;
t215 = sin(t219);
t228 = sin(qJ(1));
t248 = t228 * t215;
t216 = cos(t219);
t247 = t228 * t216;
t217 = sin(t224);
t246 = t228 * t217;
t218 = cos(t224);
t245 = t228 * t218;
t225 = sin(pkin(10));
t244 = t228 * t225;
t226 = cos(pkin(10));
t243 = t228 * t226;
t242 = t228 * t229;
t230 = cos(qJ(1));
t241 = t229 * t230;
t240 = t230 * t215;
t239 = t230 * t216;
t238 = t230 * t217;
t237 = t230 * t218;
t236 = t230 * t225;
t235 = t230 * t226;
t234 = pkin(2) * t242 + t228 * t251;
t233 = pkin(2) * t241 + t228 * pkin(7) + t230 * t251;
t232 = t227 * pkin(2) - t229 * qJ(3) + pkin(6);
t231 = g(1) * t230 + g(2) * t228;
t210 = g(1) * t228 - g(2) * t230;
t209 = g(3) * t227 + t229 * t231;
t208 = t227 * t231 - t250;
t1 = [0, -t231, t210, 0, 0, 0, 0, 0, -t209, t208, -t210, t209, -t208, -g(1) * t233 - g(2) * (-t230 * pkin(7) + t234) - g(3) * t232, -g(1) * (t227 * t236 + t243) - g(2) * (t227 * t244 - t235) + t225 * t250, -g(1) * (t227 * t235 - t244) - g(2) * (t227 * t243 + t236) + t226 * t250, -t209, -g(1) * (t228 * pkin(3) + qJ(4) * t241 + t233) - g(2) * (qJ(4) * t242 + (-pkin(3) - pkin(7)) * t230 + t234) - g(3) * (t227 * qJ(4) + t232) 0, 0, 0, 0, 0, -g(1) * (t227 * t238 + t245) - g(2) * (t227 * t246 - t237) + t217 * t250, -g(1) * (t227 * t237 - t246) - g(2) * (t227 * t245 + t238) + t218 * t250, 0, 0, 0, 0, 0, -g(1) * (t227 * t240 + t247) - g(2) * (t227 * t248 - t239) + t215 * t250, -g(1) * (t227 * t239 - t248) - g(2) * (t227 * t247 + t240) + t216 * t250;];
U_reg  = t1;
