% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:20:59
% EndTime: 2019-03-09 06:20:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (143->52), mult. (145->70), div. (0->0), fcn. (155->10), ass. (0->35)
t211 = pkin(10) + qJ(3);
t207 = sin(t211);
t233 = g(3) * t207;
t212 = qJ(4) + qJ(5);
t209 = sin(t212);
t217 = sin(qJ(1));
t232 = t217 * t209;
t210 = cos(t212);
t231 = t217 * t210;
t216 = sin(qJ(4));
t230 = t217 * t216;
t218 = cos(qJ(4));
t229 = t217 * t218;
t219 = cos(qJ(1));
t228 = t219 * t209;
t227 = t219 * t210;
t226 = t219 * t216;
t225 = t219 * t218;
t224 = pkin(4) * t216 + pkin(7) + qJ(2);
t223 = g(1) * t219 + g(2) * t217;
t206 = t218 * pkin(4) + pkin(3);
t208 = cos(t211);
t214 = cos(pkin(10));
t220 = -pkin(9) - pkin(8);
t222 = t214 * pkin(2) + t206 * t208 - t207 * t220 + pkin(1);
t199 = t208 * t232 + t227;
t201 = t208 * t228 - t231;
t221 = g(1) * t201 + g(2) * t199 + t209 * t233;
t213 = sin(pkin(10));
t203 = g(1) * t217 - g(2) * t219;
t202 = t208 * t227 + t232;
t200 = t208 * t231 - t228;
t198 = -g(3) * t208 + t207 * t223;
t197 = -g(1) * t202 - g(2) * t200 - t210 * t233;
t1 = [0, -t223, t203, -g(3) * t213 - t214 * t223, -g(3) * t214 + t213 * t223, -t203, -g(1) * (t219 * pkin(1) + t217 * qJ(2)) - g(2) * (t217 * pkin(1) - t219 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t208 * t223 - t233, t198, 0, 0, 0, 0, 0, -g(1) * (t208 * t225 + t230) - g(2) * (t208 * t229 - t226) - t218 * t233, -g(1) * (-t208 * t226 + t229) - g(2) * (-t208 * t230 - t225) + t216 * t233, 0, 0, 0, 0, 0, t197, t221, t197, -t198, -t221, -g(1) * (t202 * pkin(5) + t201 * qJ(6)) - g(2) * (t200 * pkin(5) + t199 * qJ(6)) - g(3) * (t213 * pkin(2) + t208 * t220 + pkin(6)) - (pkin(5) * t210 + qJ(6) * t209 + t206) * t233 + (-g(1) * t222 + g(2) * t224) * t219 + (-g(1) * t224 - g(2) * t222) * t217;];
U_reg  = t1;
