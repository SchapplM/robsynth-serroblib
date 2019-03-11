% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP9
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
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:27:14
% EndTime: 2019-03-09 17:27:14
% DurationCPUTime: 0.11s
% Computational Cost: add. (130->57), mult. (291->80), div. (0->0), fcn. (345->8), ass. (0->35)
t227 = sin(qJ(3));
t228 = sin(qJ(2));
t246 = t227 * t228;
t229 = sin(qJ(1));
t245 = t228 * t229;
t231 = cos(qJ(3));
t244 = t228 * t231;
t233 = cos(qJ(1));
t243 = t228 * t233;
t232 = cos(qJ(2));
t242 = t229 * t232;
t241 = t233 * t227;
t240 = t233 * t231;
t239 = t228 * pkin(2) + pkin(3) * t244 + qJ(4) * t246 + pkin(6);
t238 = g(1) * t233 + g(2) * t229;
t210 = t227 * t242 + t240;
t211 = t231 * t242 - t241;
t226 = sin(qJ(5));
t230 = cos(qJ(5));
t199 = -t210 * t230 + t211 * t226;
t212 = -t229 * t231 + t232 * t241;
t213 = t229 * t227 + t232 * t240;
t201 = -t212 * t230 + t213 * t226;
t206 = t226 * t244 - t230 * t246;
t237 = g(1) * t201 + g(2) * t199 + g(3) * t206;
t236 = t213 * pkin(3) + t229 * pkin(7) + pkin(8) * t243 + t212 * qJ(4) + (pkin(2) * t232 + pkin(1)) * t233;
t235 = g(1) * t212 + g(2) * t210 + g(3) * t246;
t234 = t229 * pkin(1) + pkin(2) * t242 + t211 * pkin(3) - t233 * pkin(7) + pkin(8) * t245 + t210 * qJ(4);
t207 = (t226 * t227 + t230 * t231) * t228;
t203 = -g(3) * t232 + t238 * t228;
t202 = t212 * t226 + t213 * t230;
t200 = t210 * t226 + t211 * t230;
t198 = -g(1) * t213 - g(2) * t211 - g(3) * t244;
t197 = -g(1) * t202 - g(2) * t200 - g(3) * t207;
t1 = [0, -t238, g(1) * t229 - g(2) * t233, 0, 0, 0, 0, 0, -g(3) * t228 - t238 * t232, t203, 0, 0, 0, 0, 0, t198, t235, t198, -t203, -t235, -g(1) * t236 - g(2) * t234 - g(3) * (-t232 * pkin(8) + t239) 0, 0, 0, 0, 0, t197, t237, t197, t203, -t237, -g(1) * (t213 * pkin(4) + t202 * pkin(5) - pkin(9) * t243 + t201 * qJ(6) + t236) - g(2) * (t211 * pkin(4) + t200 * pkin(5) - pkin(9) * t245 + t199 * qJ(6) + t234) - g(3) * (pkin(4) * t244 + t207 * pkin(5) + t206 * qJ(6) + (-pkin(8) + pkin(9)) * t232 + t239);];
U_reg  = t1;
