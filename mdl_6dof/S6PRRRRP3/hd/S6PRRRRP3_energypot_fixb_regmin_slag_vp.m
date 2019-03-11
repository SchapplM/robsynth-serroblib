% Calculate minimal parameter regressor of potential energy for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:11:06
% EndTime: 2019-03-09 00:11:06
% DurationCPUTime: 0.18s
% Computational Cost: add. (147->63), mult. (287->102), div. (0->0), fcn. (358->12), ass. (0->34)
t219 = qJ(4) + qJ(5);
t216 = sin(t219);
t224 = sin(qJ(4));
t237 = t224 * pkin(4) + pkin(5) * t216 + pkin(8);
t221 = sin(pkin(6));
t225 = sin(qJ(3));
t236 = t221 * t225;
t226 = sin(qJ(2));
t235 = t221 * t226;
t228 = cos(qJ(3));
t234 = t221 * t228;
t229 = cos(qJ(2));
t233 = t221 * t229;
t223 = cos(pkin(6));
t232 = t223 * t226;
t231 = t223 * t229;
t220 = sin(pkin(11));
t222 = cos(pkin(11));
t206 = t220 * t229 + t222 * t232;
t201 = t206 * t225 + t222 * t234;
t208 = -t220 * t232 + t222 * t229;
t203 = t208 * t225 - t220 * t234;
t209 = -t223 * t228 + t225 * t235;
t230 = g(1) * t203 + g(2) * t201 + g(3) * t209;
t227 = cos(qJ(4));
t218 = -qJ(6) - pkin(10) - pkin(9);
t217 = cos(t219);
t211 = t227 * pkin(4) + pkin(5) * t217 + pkin(3);
t210 = t223 * t225 + t226 * t234;
t207 = t220 * t231 + t222 * t226;
t205 = t220 * t226 - t222 * t231;
t204 = t208 * t228 + t220 * t236;
t202 = t206 * t228 - t222 * t236;
t1 = [-g(3) * qJ(1), 0, -g(1) * t208 - g(2) * t206 - g(3) * t235, g(1) * t207 + g(2) * t205 - g(3) * t233, 0, 0, 0, 0, 0, -g(1) * t204 - g(2) * t202 - g(3) * t210, t230, 0, 0, 0, 0, 0, -g(1) * (t204 * t227 + t207 * t224) - g(2) * (t202 * t227 + t205 * t224) - g(3) * (t210 * t227 - t224 * t233) -g(1) * (-t204 * t224 + t207 * t227) - g(2) * (-t202 * t224 + t205 * t227) - g(3) * (-t210 * t224 - t227 * t233) 0, 0, 0, 0, 0, -g(1) * (t204 * t217 + t207 * t216) - g(2) * (t202 * t217 + t205 * t216) - g(3) * (t210 * t217 - t216 * t233) -g(1) * (-t204 * t216 + t207 * t217) - g(2) * (-t202 * t216 + t205 * t217) - g(3) * (-t210 * t216 - t217 * t233) -t230, -g(1) * (t222 * pkin(1) + t208 * pkin(2) - t203 * t218 + t204 * t211 + t237 * t207) - g(2) * (t220 * pkin(1) + t206 * pkin(2) - t201 * t218 + t202 * t211 + t237 * t205) - g(3) * (t223 * pkin(7) - t209 * t218 + t210 * t211 + qJ(1)) + (-g(3) * (pkin(2) * t226 - t237 * t229) + (-g(1) * t220 + g(2) * t222) * pkin(7)) * t221;];
U_reg  = t1;
