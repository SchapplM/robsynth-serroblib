% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:51:32
% EndTime: 2019-03-09 16:51:33
% DurationCPUTime: 0.14s
% Computational Cost: add. (153->59), mult. (157->75), div. (0->0), fcn. (167->10), ass. (0->35)
t225 = sin(qJ(2));
t241 = g(3) * t225;
t224 = sin(qJ(3));
t240 = t224 * pkin(3);
t223 = -qJ(4) - pkin(8);
t227 = cos(qJ(3));
t214 = pkin(3) * t227 + pkin(2);
t239 = t223 * t225;
t226 = sin(qJ(1));
t238 = t226 * t224;
t228 = cos(qJ(2));
t237 = t226 * t228;
t222 = qJ(3) + pkin(10);
t215 = qJ(5) + t222;
t212 = sin(t215);
t229 = cos(qJ(1));
t236 = t229 * t212;
t213 = cos(t215);
t235 = t229 * t213;
t234 = t229 * t224;
t233 = t229 * t227;
t232 = pkin(1) * t229 + pkin(7) * t226;
t231 = g(1) * t229 + g(2) * t226;
t204 = t212 * t237 + t235;
t206 = -t213 * t226 + t228 * t236;
t230 = g(1) * t206 + g(2) * t204 + t212 * t241;
t221 = -pkin(9) + t223;
t217 = t226 * pkin(1);
t210 = pkin(4) * sin(t222) + t240;
t209 = pkin(4) * cos(t222) + t214;
t208 = -g(3) * t228 + t225 * t231;
t207 = t212 * t226 + t228 * t235;
t205 = t213 * t237 - t236;
t203 = -g(1) * t207 - g(2) * t205 - t213 * t241;
t1 = [0, -t231, g(1) * t226 - g(2) * t229, 0, 0, 0, 0, 0, -t228 * t231 - t241, t208, 0, 0, 0, 0, 0, -g(1) * (t228 * t233 + t238) - g(2) * (t227 * t237 - t234) - t227 * t241, -g(1) * (t226 * t227 - t228 * t234) - g(2) * (-t224 * t237 - t233) + t224 * t241, -t208, -g(1) * (pkin(3) * t238 + t232) - g(2) * (t214 * t237 - t226 * t239 + t217) - g(3) * (t214 * t225 + t223 * t228 + pkin(6)) + (-g(1) * (t214 * t228 - t239) - g(2) * (-pkin(7) - t240)) * t229, 0, 0, 0, 0, 0, t203, t230, t203, -t208, -t230, -g(1) * (t229 * t228 * t209 + t207 * pkin(5) + t206 * qJ(6) + t226 * t210 + t232) - g(2) * (t205 * pkin(5) + t204 * qJ(6) + t209 * t237 + t217 + (-pkin(7) - t210) * t229) - g(3) * (t228 * t221 + pkin(6)) + (-g(3) * (pkin(5) * t213 + qJ(6) * t212 + t209) + t231 * t221) * t225;];
U_reg  = t1;
