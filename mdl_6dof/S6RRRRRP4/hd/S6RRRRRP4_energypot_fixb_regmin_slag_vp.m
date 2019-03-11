% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP4
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
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:15:40
% EndTime: 2019-03-10 01:15:40
% DurationCPUTime: 0.08s
% Computational Cost: add. (138->47), mult. (136->63), div. (0->0), fcn. (149->10), ass. (0->34)
t213 = qJ(2) + qJ(3);
t209 = sin(t213);
t234 = g(3) * t209;
t212 = qJ(4) + qJ(5);
t208 = sin(t212);
t216 = sin(qJ(1));
t233 = t216 * t208;
t210 = cos(t212);
t232 = t216 * t210;
t214 = sin(qJ(4));
t231 = t216 * t214;
t217 = cos(qJ(4));
t230 = t216 * t217;
t219 = cos(qJ(1));
t229 = t219 * t208;
t228 = t219 * t210;
t227 = t219 * t214;
t226 = t219 * t217;
t225 = pkin(4) * t214 + pkin(7) + pkin(8);
t224 = g(1) * t219 + g(2) * t216;
t206 = t217 * pkin(4) + pkin(3);
t211 = cos(t213);
t218 = cos(qJ(2));
t220 = -pkin(10) - pkin(9);
t223 = t218 * pkin(2) + t206 * t211 - t209 * t220 + pkin(1);
t201 = t211 * t233 + t228;
t203 = t211 * t229 - t232;
t222 = g(1) * t203 + g(2) * t201 + t208 * t234;
t215 = sin(qJ(2));
t204 = t211 * t228 + t233;
t202 = t211 * t232 - t229;
t200 = -g(3) * t211 + t224 * t209;
t199 = -g(1) * t204 - g(2) * t202 - t210 * t234;
t1 = [0, -t224, g(1) * t216 - g(2) * t219, 0, 0, 0, 0, 0, -g(3) * t215 - t224 * t218, -g(3) * t218 + t224 * t215, 0, 0, 0, 0, 0, -t224 * t211 - t234, t200, 0, 0, 0, 0, 0, -g(1) * (t211 * t226 + t231) - g(2) * (t211 * t230 - t227) - t217 * t234, -g(1) * (-t211 * t227 + t230) - g(2) * (-t211 * t231 - t226) + t214 * t234, 0, 0, 0, 0, 0, t199, t222, t199, -t200, -t222, -g(1) * (t204 * pkin(5) + t203 * qJ(6)) - g(2) * (t202 * pkin(5) + t201 * qJ(6)) - g(3) * (t215 * pkin(2) + t211 * t220 + pkin(6)) - (pkin(5) * t210 + qJ(6) * t208 + t206) * t234 + (-g(1) * t223 + g(2) * t225) * t219 + (-g(1) * t225 - g(2) * t223) * t216;];
U_reg  = t1;
