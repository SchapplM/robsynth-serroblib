% Calculate minimal parameter regressor of potential energy for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:39:50
% EndTime: 2019-03-09 00:39:50
% DurationCPUTime: 0.13s
% Computational Cost: add. (135->51), mult. (191->95), div. (0->0), fcn. (244->14), ass. (0->32)
t224 = sin(pkin(12));
t225 = sin(pkin(6));
t241 = t224 * t225;
t226 = cos(pkin(12));
t240 = t225 * t226;
t229 = sin(qJ(3));
t239 = t225 * t229;
t230 = sin(qJ(2));
t238 = t225 * t230;
t232 = cos(qJ(3));
t237 = t225 * t232;
t233 = cos(qJ(2));
t236 = t225 * t233;
t227 = cos(pkin(6));
t235 = t227 * t230;
t234 = t227 * t233;
t223 = qJ(3) + qJ(4);
t231 = cos(qJ(6));
t228 = sin(qJ(6));
t222 = qJ(5) + t223;
t221 = cos(t223);
t220 = sin(t223);
t219 = cos(t222);
t218 = sin(t222);
t217 = -t224 * t235 + t226 * t233;
t216 = t224 * t234 + t226 * t230;
t215 = t224 * t233 + t226 * t235;
t214 = t224 * t230 - t226 * t234;
t213 = t227 * t218 + t219 * t238;
t212 = t217 * t219 + t218 * t241;
t211 = t215 * t219 - t218 * t240;
t1 = [-g(3) * qJ(1), 0, -g(1) * t217 - g(2) * t215 - g(3) * t238, g(1) * t216 + g(2) * t214 - g(3) * t236, 0, 0, 0, 0, 0, -g(1) * (t217 * t232 + t224 * t239) - g(2) * (t215 * t232 - t226 * t239) - g(3) * (t227 * t229 + t230 * t237) -g(1) * (-t217 * t229 + t224 * t237) - g(2) * (-t215 * t229 - t226 * t237) - g(3) * (t227 * t232 - t229 * t238) 0, 0, 0, 0, 0, -g(1) * (t217 * t221 + t220 * t241) - g(2) * (t215 * t221 - t220 * t240) - g(3) * (t227 * t220 + t221 * t238) -g(1) * (-t217 * t220 + t221 * t241) - g(2) * (-t215 * t220 - t221 * t240) - g(3) * (-t220 * t238 + t227 * t221) 0, 0, 0, 0, 0, -g(1) * t212 - g(2) * t211 - g(3) * t213, -g(1) * (-t217 * t218 + t219 * t241) - g(2) * (-t215 * t218 - t219 * t240) - g(3) * (-t218 * t238 + t227 * t219) 0, 0, 0, 0, 0, -g(1) * (t212 * t231 + t216 * t228) - g(2) * (t211 * t231 + t214 * t228) - g(3) * (t213 * t231 - t228 * t236) -g(1) * (-t212 * t228 + t216 * t231) - g(2) * (-t211 * t228 + t214 * t231) - g(3) * (-t213 * t228 - t231 * t236);];
U_reg  = t1;
