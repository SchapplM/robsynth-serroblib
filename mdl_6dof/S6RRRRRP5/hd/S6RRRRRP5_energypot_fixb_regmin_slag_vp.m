% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:23:03
% EndTime: 2019-03-10 01:23:03
% DurationCPUTime: 0.12s
% Computational Cost: add. (111->47), mult. (113->68), div. (0->0), fcn. (122->10), ass. (0->27)
t214 = sin(qJ(2));
t229 = g(3) * t214;
t212 = qJ(3) + qJ(4);
t211 = qJ(5) + t212;
t206 = sin(t211);
t208 = sin(t212);
t213 = sin(qJ(3));
t228 = t213 * pkin(3) + pkin(4) * t208 + pkin(5) * t206 + pkin(7);
t215 = sin(qJ(1));
t217 = cos(qJ(2));
t227 = t215 * t217;
t218 = cos(qJ(1));
t226 = t218 * t206;
t207 = cos(t211);
t225 = t218 * t207;
t224 = t218 * t208;
t209 = cos(t212);
t223 = t218 * t209;
t222 = t218 * t213;
t216 = cos(qJ(3));
t221 = t218 * t216;
t220 = g(1) * t218 + g(2) * t215;
t203 = t216 * pkin(3) + pkin(4) * t209 + pkin(5) * t207 + pkin(2);
t210 = -qJ(6) - pkin(10) - pkin(9) - pkin(8);
t219 = t203 * t217 - t210 * t214 + pkin(1);
t204 = -g(3) * t217 + t220 * t214;
t1 = [0, -t220, g(1) * t215 - g(2) * t218, 0, 0, 0, 0, 0, -t220 * t217 - t229, t204, 0, 0, 0, 0, 0, -g(1) * (t215 * t213 + t217 * t221) - g(2) * (t216 * t227 - t222) - t216 * t229, -g(1) * (t215 * t216 - t217 * t222) - g(2) * (-t213 * t227 - t221) + t213 * t229, 0, 0, 0, 0, 0, -g(1) * (t215 * t208 + t217 * t223) - g(2) * (t209 * t227 - t224) - t209 * t229, -g(1) * (t215 * t209 - t217 * t224) - g(2) * (-t208 * t227 - t223) + t208 * t229, 0, 0, 0, 0, 0, -g(1) * (t215 * t206 + t217 * t225) - g(2) * (t207 * t227 - t226) - t207 * t229, -g(1) * (t215 * t207 - t217 * t226) - g(2) * (-t206 * t227 - t225) + t206 * t229, -t204, -g(3) * (t214 * t203 + t217 * t210 + pkin(6)) + (-g(1) * t219 + g(2) * t228) * t218 + (-g(1) * t228 - g(2) * t219) * t215;];
U_reg  = t1;
