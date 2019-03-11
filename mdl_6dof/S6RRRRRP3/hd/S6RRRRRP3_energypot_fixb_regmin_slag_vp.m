% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP3
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
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:09:29
% EndTime: 2019-03-10 01:09:30
% DurationCPUTime: 0.08s
% Computational Cost: add. (101->42), mult. (101->57), div. (0->0), fcn. (106->10), ass. (0->28)
t206 = qJ(2) + qJ(3);
t201 = sin(t206);
t225 = g(3) * t201;
t205 = qJ(4) + qJ(5);
t200 = sin(t205);
t209 = sin(qJ(1));
t224 = t209 * t200;
t202 = cos(t205);
t223 = t209 * t202;
t207 = sin(qJ(4));
t222 = t209 * t207;
t210 = cos(qJ(4));
t221 = t209 * t210;
t212 = cos(qJ(1));
t220 = t212 * t200;
t219 = t212 * t202;
t218 = t212 * t207;
t217 = t212 * t210;
t216 = t207 * pkin(4) + pkin(5) * t200 + pkin(7) + pkin(8);
t215 = g(1) * t212 + g(2) * t209;
t197 = t210 * pkin(4) + pkin(5) * t202 + pkin(3);
t203 = cos(t206);
t204 = -qJ(6) - pkin(10) - pkin(9);
t211 = cos(qJ(2));
t214 = t211 * pkin(2) + t197 * t203 - t201 * t204 + pkin(1);
t208 = sin(qJ(2));
t196 = -g(3) * t203 + t215 * t201;
t1 = [0, -t215, g(1) * t209 - g(2) * t212, 0, 0, 0, 0, 0, -g(3) * t208 - t215 * t211, -g(3) * t211 + t215 * t208, 0, 0, 0, 0, 0, -t215 * t203 - t225, t196, 0, 0, 0, 0, 0, -g(1) * (t203 * t217 + t222) - g(2) * (t203 * t221 - t218) - t210 * t225, -g(1) * (-t203 * t218 + t221) - g(2) * (-t203 * t222 - t217) + t207 * t225, 0, 0, 0, 0, 0, -g(1) * (t203 * t219 + t224) - g(2) * (t203 * t223 - t220) - t202 * t225, -g(1) * (-t203 * t220 + t223) - g(2) * (-t203 * t224 - t219) + t200 * t225, -t196, -g(3) * (t208 * pkin(2) + t201 * t197 + t203 * t204 + pkin(6)) + (-g(1) * t214 + g(2) * t216) * t212 + (-g(1) * t216 - g(2) * t214) * t209;];
U_reg  = t1;
