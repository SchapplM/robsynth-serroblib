% Calculate minimal parameter regressor of potential energy for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:42
% EndTime: 2019-03-09 15:26:42
% DurationCPUTime: 0.07s
% Computational Cost: add. (106->37), mult. (94->50), div. (0->0), fcn. (92->10), ass. (0->26)
t209 = qJ(2) + qJ(3);
t203 = pkin(10) + t209;
t201 = cos(t203);
t225 = g(3) * t201;
t210 = sin(qJ(6));
t212 = sin(qJ(1));
t224 = t212 * t210;
t213 = cos(qJ(6));
t223 = t212 * t213;
t215 = cos(qJ(1));
t222 = t215 * t210;
t221 = t215 * t213;
t205 = cos(t209);
t214 = cos(qJ(2));
t197 = t214 * pkin(2) + pkin(3) * t205 + pkin(1);
t208 = -qJ(4) - pkin(8) - pkin(7);
t220 = t212 * t197 + t215 * t208;
t204 = sin(t209);
t211 = sin(qJ(2));
t219 = t211 * pkin(2) + pkin(3) * t204 + pkin(6);
t218 = t215 * t197 - t212 * t208;
t217 = g(1) * t215 + g(2) * t212;
t200 = sin(t203);
t216 = pkin(4) * t201 + qJ(5) * t200;
t198 = g(1) * t212 - g(2) * t215;
t1 = [0, -t217, t198, 0, 0, 0, 0, 0, -g(3) * t211 - t217 * t214, -g(3) * t214 + t217 * t211, 0, 0, 0, 0, 0, -g(3) * t204 - t217 * t205, -g(3) * t205 + t217 * t204, -t198, -g(1) * t218 - g(2) * t220 - g(3) * t219, -t198, g(3) * t200 + t217 * t201, -t217 * t200 + t225, -g(1) * (t216 * t215 + t218) - g(2) * (t216 * t212 + t220) - g(3) * (t200 * pkin(4) - t201 * qJ(5) + t219) 0, 0, 0, 0, 0, -g(1) * (t200 * t222 + t223) - g(2) * (t200 * t224 - t221) + t210 * t225, -g(1) * (t200 * t221 - t224) - g(2) * (t200 * t223 + t222) + t213 * t225;];
U_reg  = t1;
