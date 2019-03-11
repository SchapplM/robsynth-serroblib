% Calculate minimal parameter regressor of potential energy for
% S6PRRRRR2
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
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:45:22
% EndTime: 2019-03-09 00:45:22
% DurationCPUTime: 0.10s
% Computational Cost: add. (133->51), mult. (217->95), div. (0->0), fcn. (280->14), ass. (0->32)
t227 = sin(pkin(12));
t228 = sin(pkin(6));
t244 = t227 * t228;
t229 = cos(pkin(12));
t243 = t228 * t229;
t232 = sin(qJ(3));
t242 = t228 * t232;
t233 = sin(qJ(2));
t241 = t228 * t233;
t235 = cos(qJ(3));
t240 = t228 * t235;
t236 = cos(qJ(2));
t239 = t228 * t236;
t230 = cos(pkin(6));
t238 = t230 * t233;
t237 = t230 * t236;
t234 = cos(qJ(5));
t231 = sin(qJ(5));
t226 = qJ(3) + qJ(4);
t225 = qJ(5) + qJ(6);
t224 = cos(t226);
t223 = cos(t225);
t222 = sin(t226);
t221 = sin(t225);
t220 = -t227 * t238 + t229 * t236;
t219 = t227 * t237 + t229 * t233;
t218 = t227 * t236 + t229 * t238;
t217 = t227 * t233 - t229 * t237;
t216 = t230 * t222 + t224 * t241;
t215 = t220 * t224 + t222 * t244;
t214 = t218 * t224 - t222 * t243;
t1 = [-g(3) * qJ(1), 0, -g(1) * t220 - g(2) * t218 - g(3) * t241, g(1) * t219 + g(2) * t217 - g(3) * t239, 0, 0, 0, 0, 0, -g(1) * (t220 * t235 + t227 * t242) - g(2) * (t218 * t235 - t229 * t242) - g(3) * (t230 * t232 + t233 * t240) -g(1) * (-t220 * t232 + t227 * t240) - g(2) * (-t218 * t232 - t229 * t240) - g(3) * (t230 * t235 - t232 * t241) 0, 0, 0, 0, 0, -g(1) * t215 - g(2) * t214 - g(3) * t216, -g(1) * (-t220 * t222 + t224 * t244) - g(2) * (-t218 * t222 - t224 * t243) - g(3) * (-t222 * t241 + t230 * t224) 0, 0, 0, 0, 0, -g(1) * (t215 * t234 + t219 * t231) - g(2) * (t214 * t234 + t217 * t231) - g(3) * (t216 * t234 - t231 * t239) -g(1) * (-t215 * t231 + t219 * t234) - g(2) * (-t214 * t231 + t217 * t234) - g(3) * (-t216 * t231 - t234 * t239) 0, 0, 0, 0, 0, -g(1) * (t215 * t223 + t219 * t221) - g(2) * (t214 * t223 + t217 * t221) - g(3) * (t216 * t223 - t221 * t239) -g(1) * (-t215 * t221 + t219 * t223) - g(2) * (-t214 * t221 + t217 * t223) - g(3) * (-t216 * t221 - t223 * t239);];
U_reg  = t1;
