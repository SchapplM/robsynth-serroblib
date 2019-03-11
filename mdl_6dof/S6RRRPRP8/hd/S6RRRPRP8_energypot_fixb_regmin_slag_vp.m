% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:19:39
% EndTime: 2019-03-09 17:19:39
% DurationCPUTime: 0.11s
% Computational Cost: add. (103->56), mult. (215->79), div. (0->0), fcn. (243->8), ass. (0->32)
t214 = sin(qJ(2));
t234 = g(3) * t214;
t212 = sin(qJ(5));
t213 = sin(qJ(3));
t233 = t212 * t213;
t232 = t213 * t214;
t215 = sin(qJ(1));
t231 = t214 * t215;
t217 = cos(qJ(3));
t230 = t214 * t217;
t219 = cos(qJ(1));
t229 = t214 * t219;
t218 = cos(qJ(2));
t228 = t215 * t218;
t227 = t219 * t213;
t226 = t219 * t217;
t225 = pkin(5) * t212 + qJ(4);
t224 = t214 * pkin(2) + pkin(3) * t230 + qJ(4) * t232 + pkin(6);
t223 = g(1) * t219 + g(2) * t215;
t198 = t215 * t213 + t218 * t226;
t222 = t198 * pkin(3) + t215 * pkin(7) + pkin(8) * t229 + (pkin(2) * t218 + pkin(1)) * t219;
t196 = t217 * t228 - t227;
t221 = t215 * pkin(1) + pkin(2) * t228 + t196 * pkin(3) - t219 * pkin(7) + pkin(8) * t231;
t195 = t213 * t228 + t226;
t197 = -t215 * t217 + t218 * t227;
t220 = g(1) * t197 + g(2) * t195 + g(3) * t232;
t216 = cos(qJ(5));
t211 = -qJ(6) - pkin(9);
t206 = t216 * pkin(5) + pkin(4);
t192 = -g(3) * t218 + t223 * t214;
t191 = -g(1) * t198 - g(2) * t196 - g(3) * t230;
t1 = [0, -t223, g(1) * t215 - g(2) * t219, 0, 0, 0, 0, 0, -t223 * t218 - t234, t192, 0, 0, 0, 0, 0, t191, t220, t191, -t192, -t220, -g(1) * (t197 * qJ(4) + t222) - g(2) * (t195 * qJ(4) + t221) - g(3) * (-t218 * pkin(8) + t224) 0, 0, 0, 0, 0, -g(1) * (t197 * t212 + t198 * t216) - g(2) * (t195 * t212 + t196 * t216) - (t216 * t217 + t233) * t234, -g(1) * (t197 * t216 - t198 * t212) - g(2) * (t195 * t216 - t196 * t212) - (-t212 * t217 + t213 * t216) * t234, t192, -g(1) * (t225 * t197 + t198 * t206 + t211 * t229 + t222) - g(2) * (t225 * t195 + t196 * t206 + t211 * t231 + t221) - g(3) * ((-pkin(8) - t211) * t218 + (pkin(5) * t233 + t206 * t217) * t214 + t224);];
U_reg  = t1;
