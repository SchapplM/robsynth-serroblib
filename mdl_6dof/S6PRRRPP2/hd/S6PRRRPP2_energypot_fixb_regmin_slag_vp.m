% Calculate minimal parameter regressor of potential energy for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRPP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:52:40
% EndTime: 2019-03-08 22:52:40
% DurationCPUTime: 0.12s
% Computational Cost: add. (244->64), mult. (584->93), div. (0->0), fcn. (744->10), ass. (0->41)
t227 = sin(pkin(6));
t247 = pkin(7) * t227;
t246 = pkin(9) - qJ(6);
t231 = sin(qJ(3));
t245 = t227 * t231;
t232 = sin(qJ(2));
t244 = t227 * t232;
t234 = cos(qJ(3));
t243 = t227 * t234;
t235 = cos(qJ(2));
t242 = t227 * t235;
t229 = cos(pkin(6));
t241 = t229 * t232;
t240 = t229 * t235;
t226 = sin(pkin(10));
t228 = cos(pkin(10));
t213 = t226 * t235 + t228 * t241;
t202 = t213 * t234 - t228 * t245;
t212 = t226 * t232 - t228 * t240;
t230 = sin(qJ(4));
t233 = cos(qJ(4));
t194 = t202 * t230 - t212 * t233;
t215 = -t226 * t241 + t228 * t235;
t204 = t215 * t234 + t226 * t245;
t214 = t226 * t240 + t228 * t232;
t196 = t204 * t230 - t214 * t233;
t217 = t229 * t231 + t232 * t243;
t205 = t217 * t230 + t233 * t242;
t239 = g(1) * t196 + g(2) * t194 + g(3) * t205;
t201 = t213 * t231 + t228 * t243;
t203 = t215 * t231 - t226 * t243;
t216 = -t229 * t234 + t231 * t244;
t191 = g(1) * t203 + g(2) * t201 + g(3) * t216;
t197 = t204 * t233 + t214 * t230;
t238 = t228 * pkin(1) + t215 * pkin(2) + t204 * pkin(3) + t197 * pkin(4) + t214 * pkin(8) + t196 * qJ(5) + t226 * t247;
t206 = t217 * t233 - t230 * t242;
t237 = pkin(2) * t244 + t217 * pkin(3) + t206 * pkin(4) + t229 * pkin(7) - pkin(8) * t242 + t205 * qJ(5) + qJ(1);
t195 = t202 * t233 + t212 * t230;
t236 = t226 * pkin(1) + t213 * pkin(2) + t202 * pkin(3) + t195 * pkin(4) + t212 * pkin(8) + t194 * qJ(5) - t228 * t247;
t190 = -g(1) * t197 - g(2) * t195 - g(3) * t206;
t1 = [-g(3) * qJ(1), 0, -g(1) * t215 - g(2) * t213 - g(3) * t244, g(1) * t214 + g(2) * t212 - g(3) * t242, 0, 0, 0, 0, 0, -g(1) * t204 - g(2) * t202 - g(3) * t217, t191, 0, 0, 0, 0, 0, t190, t239, t190, -t191, -t239, -g(1) * (t203 * pkin(9) + t238) - g(2) * (t201 * pkin(9) + t236) - g(3) * (t216 * pkin(9) + t237) t190, -t239, t191, -g(1) * (t197 * pkin(5) + t246 * t203 + t238) - g(2) * (t195 * pkin(5) + t246 * t201 + t236) - g(3) * (t206 * pkin(5) + t246 * t216 + t237);];
U_reg  = t1;
