% Calculate minimal parameter regressor of potential energy for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:33:00
% EndTime: 2019-03-09 16:33:00
% DurationCPUTime: 0.09s
% Computational Cost: add. (105->40), mult. (94->54), div. (0->0), fcn. (92->10), ass. (0->28)
t209 = qJ(2) + qJ(3);
t204 = pkin(10) + t209;
t200 = sin(t204);
t201 = cos(t204);
t214 = cos(qJ(5));
t203 = t214 * pkin(5) + pkin(4);
t210 = -qJ(6) - pkin(9);
t227 = -t200 * t210 + t201 * t203;
t226 = g(3) * t200;
t211 = sin(qJ(5));
t213 = sin(qJ(1));
t223 = t213 * t211;
t222 = t213 * t214;
t216 = cos(qJ(1));
t221 = t216 * t211;
t220 = t216 * t214;
t206 = cos(t209);
t215 = cos(qJ(2));
t197 = t215 * pkin(2) + pkin(3) * t206 + pkin(1);
t208 = -qJ(4) - pkin(8) - pkin(7);
t219 = t213 * t197 + t216 * t208;
t205 = sin(t209);
t212 = sin(qJ(2));
t218 = t212 * pkin(2) + pkin(3) * t205 + pkin(6);
t217 = g(1) * t216 + g(2) * t213;
t198 = g(1) * t213 - g(2) * t216;
t196 = t216 * t197;
t1 = [0, -t217, t198, 0, 0, 0, 0, 0, -g(3) * t212 - t217 * t215, -g(3) * t215 + t217 * t212, 0, 0, 0, 0, 0, -g(3) * t205 - t217 * t206, -g(3) * t206 + t217 * t205, -t198, -g(1) * (-t213 * t208 + t196) - g(2) * t219 - g(3) * t218, 0, 0, 0, 0, 0, -g(1) * (t201 * t220 + t223) - g(2) * (t201 * t222 - t221) - t214 * t226, -g(1) * (-t201 * t221 + t222) - g(2) * (-t201 * t223 - t220) + t211 * t226, g(3) * t201 - t217 * t200, -g(1) * (t227 * t216 + t196) - g(2) * (-pkin(5) * t221 + t219) - g(3) * (t200 * t203 + t201 * t210 + t218) + (-g(1) * (pkin(5) * t211 - t208) - g(2) * t227) * t213;];
U_reg  = t1;
