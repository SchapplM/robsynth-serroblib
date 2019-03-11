% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:13:22
% EndTime: 2019-03-09 18:13:22
% DurationCPUTime: 0.08s
% Computational Cost: add. (96->34), mult. (118->50), div. (0->0), fcn. (131->10), ass. (0->23)
t216 = sin(qJ(1));
t220 = cos(qJ(1));
t224 = g(1) * t220 + g(2) * t216;
t212 = qJ(2) + qJ(3);
t210 = sin(t212);
t211 = cos(t212);
t214 = sin(qJ(5));
t218 = cos(qJ(5));
t207 = t210 * t218 - t211 * t214;
t225 = g(3) * t207;
t223 = t210 * t214 + t211 * t218;
t219 = cos(qJ(2));
t222 = t219 * pkin(2) + pkin(3) * t211 + qJ(4) * t210 + pkin(1);
t221 = -pkin(8) - pkin(7);
t217 = cos(qJ(6));
t215 = sin(qJ(2));
t213 = sin(qJ(6));
t208 = g(1) * t216 - g(2) * t220;
t206 = t223 * t220;
t205 = t223 * t216;
t204 = -g(3) * t210 - t224 * t211;
t203 = -g(3) * t211 + t224 * t210;
t1 = [0, -t224, t208, 0, 0, 0, 0, 0, -g(3) * t215 - t224 * t219, -g(3) * t219 + t224 * t215, 0, 0, 0, 0, 0, t204, t203, t204, -t208, -t203, -g(3) * (t215 * pkin(2) + t210 * pkin(3) - t211 * qJ(4) + pkin(6)) + (-g(1) * t222 - g(2) * t221) * t220 + (g(1) * t221 - g(2) * t222) * t216, 0, 0, 0, 0, 0, -g(1) * t206 - g(2) * t205 - t225, g(3) * t223 - t224 * t207, 0, 0, 0, 0, 0, -g(1) * (t206 * t217 - t216 * t213) - g(2) * (t205 * t217 + t220 * t213) - t217 * t225, -g(1) * (-t206 * t213 - t216 * t217) - g(2) * (-t205 * t213 + t220 * t217) + t213 * t225;];
U_reg  = t1;
