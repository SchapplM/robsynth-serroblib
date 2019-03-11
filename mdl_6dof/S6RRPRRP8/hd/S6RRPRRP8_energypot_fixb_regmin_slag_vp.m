% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:25:15
% EndTime: 2019-03-09 12:25:15
% DurationCPUTime: 0.14s
% Computational Cost: add. (163->62), mult. (170->86), div. (0->0), fcn. (184->10), ass. (0->35)
t225 = sin(qJ(2));
t241 = g(3) * t225;
t226 = sin(qJ(1));
t227 = cos(qJ(2));
t240 = t226 * t227;
t222 = pkin(10) + qJ(4);
t216 = qJ(5) + t222;
t212 = sin(t216);
t228 = cos(qJ(1));
t239 = t228 * t212;
t213 = cos(t216);
t238 = t228 * t213;
t214 = sin(t222);
t237 = t228 * t214;
t215 = cos(t222);
t236 = t228 * t215;
t223 = sin(pkin(10));
t235 = t228 * t223;
t224 = cos(pkin(10));
t234 = t228 * t224;
t233 = t228 * pkin(1) + t226 * pkin(7);
t232 = t226 * pkin(1) - t228 * pkin(7);
t231 = g(1) * t228 + g(2) * t226;
t230 = pkin(2) * t227 + qJ(3) * t225;
t204 = t212 * t240 + t238;
t206 = -t226 * t213 + t227 * t239;
t229 = g(1) * t206 + g(2) * t204 + t212 * t241;
t221 = -pkin(9) - pkin(8) - qJ(3);
t210 = t223 * pkin(3) + pkin(4) * t214;
t209 = t224 * pkin(3) + pkin(4) * t215 + pkin(2);
t208 = -g(3) * t227 + t231 * t225;
t207 = t226 * t212 + t227 * t238;
t205 = t213 * t240 - t239;
t203 = -g(1) * t207 - g(2) * t205 - t213 * t241;
t1 = [0, -t231, g(1) * t226 - g(2) * t228, 0, 0, 0, 0, 0, -t231 * t227 - t241, t208, -g(1) * (t226 * t223 + t227 * t234) - g(2) * (t224 * t240 - t235) - t224 * t241, -g(1) * (t226 * t224 - t227 * t235) - g(2) * (-t223 * t240 - t234) + t223 * t241, -t208, -g(1) * (t230 * t228 + t233) - g(2) * (t230 * t226 + t232) - g(3) * (t225 * pkin(2) - t227 * qJ(3) + pkin(6)) 0, 0, 0, 0, 0, -g(1) * (t226 * t214 + t227 * t236) - g(2) * (t215 * t240 - t237) - t215 * t241, -g(1) * (t226 * t215 - t227 * t237) - g(2) * (-t214 * t240 - t236) + t214 * t241, 0, 0, 0, 0, 0, t203, t229, t203, -t208, -t229, -g(1) * (t228 * t227 * t209 + t207 * pkin(5) + t206 * qJ(6) + t226 * t210 + t233) - g(2) * (t205 * pkin(5) + t204 * qJ(6) + t209 * t240 - t228 * t210 + t232) - g(3) * (t227 * t221 + pkin(6)) + (-g(3) * (pkin(5) * t213 + qJ(6) * t212 + t209) + t231 * t221) * t225;];
U_reg  = t1;
