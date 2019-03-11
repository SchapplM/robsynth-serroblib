% Calculate minimal parameter regressor of potential energy for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:22
% EndTime: 2019-03-09 02:48:22
% DurationCPUTime: 0.12s
% Computational Cost: add. (147->53), mult. (189->74), div. (0->0), fcn. (206->10), ass. (0->31)
t226 = pkin(9) + qJ(3);
t223 = sin(t226);
t224 = cos(t226);
t230 = cos(pkin(9));
t248 = t230 * pkin(2) + pkin(3) * t224 + qJ(4) * t223 + pkin(1);
t246 = g(3) * t223;
t227 = sin(pkin(10));
t233 = sin(qJ(1));
t244 = t233 * t227;
t229 = cos(pkin(10));
t243 = t233 * t229;
t235 = cos(qJ(1));
t242 = t235 * t227;
t241 = t235 * t229;
t231 = -pkin(7) - qJ(2);
t240 = t235 * t231 + t248 * t233;
t239 = g(1) * t235 + g(2) * t233;
t228 = sin(pkin(9));
t238 = t228 * pkin(2) + t223 * pkin(3) - t224 * qJ(4) + pkin(6);
t237 = -t233 * t231 + t248 * t235;
t207 = t224 * t244 + t241;
t209 = t224 * t242 - t243;
t236 = g(1) * t209 + g(2) * t207 + t227 * t246;
t234 = cos(qJ(6));
t232 = sin(qJ(6));
t215 = g(1) * t233 - g(2) * t235;
t210 = t224 * t241 + t244;
t208 = t224 * t243 - t242;
t206 = -g(3) * t224 + t239 * t223;
t205 = -g(1) * t210 - g(2) * t208 - t229 * t246;
t1 = [0, -t239, t215, -g(3) * t228 - t239 * t230, -g(3) * t230 + t239 * t228, -t215, -g(1) * (t235 * pkin(1) + t233 * qJ(2)) - g(2) * (t233 * pkin(1) - t235 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t239 * t224 - t246, t206, t205, t236, -t206, -g(1) * t237 - g(2) * t240 - g(3) * t238, t205, -t206, -t236, -g(1) * (t210 * pkin(4) + t209 * qJ(5) + t237) - g(2) * (t208 * pkin(4) + t207 * qJ(5) + t240) - g(3) * ((pkin(4) * t229 + qJ(5) * t227) * t223 + t238) 0, 0, 0, 0, 0, -g(1) * (t209 * t232 + t210 * t234) - g(2) * (t207 * t232 + t208 * t234) - (t227 * t232 + t229 * t234) * t246, -g(1) * (t209 * t234 - t210 * t232) - g(2) * (t207 * t234 - t208 * t232) - (t227 * t234 - t229 * t232) * t246;];
U_reg  = t1;
