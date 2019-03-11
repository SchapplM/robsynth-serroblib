% Calculate minimal parameter regressor of potential energy for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:46
% EndTime: 2019-03-08 20:48:46
% DurationCPUTime: 0.12s
% Computational Cost: add. (106->53), mult. (231->91), div. (0->0), fcn. (289->12), ass. (0->30)
t207 = sin(pkin(6));
t224 = pkin(7) * t207;
t211 = sin(qJ(4));
t223 = t207 * t211;
t212 = sin(qJ(2));
t222 = t207 * t212;
t214 = cos(qJ(4));
t221 = t207 * t214;
t215 = cos(qJ(2));
t220 = t207 * t215;
t209 = cos(pkin(6));
t219 = t209 * t212;
t218 = t209 * t215;
t206 = sin(pkin(11));
t208 = cos(pkin(11));
t197 = t206 * t212 - t208 * t218;
t199 = t206 * t218 + t208 * t212;
t217 = -g(1) * t199 - g(2) * t197 + g(3) * t220;
t198 = t206 * t215 + t208 * t219;
t200 = -t206 * t219 + t208 * t215;
t216 = g(1) * t200 + g(2) * t198 + g(3) * t222;
t213 = cos(qJ(5));
t210 = sin(qJ(5));
t205 = qJ(5) + qJ(6);
t204 = cos(t205);
t203 = sin(t205);
t201 = t209 * t214 - t211 * t220;
t196 = t197 * t211 - t208 * t221;
t195 = t199 * t211 + t206 * t221;
t1 = [-g(3) * qJ(1), 0, -t216, -t217, t216, t217, -g(1) * (t208 * pkin(1) + t200 * pkin(2) + t199 * qJ(3) + t206 * t224) - g(2) * (t206 * pkin(1) + t198 * pkin(2) + t197 * qJ(3) - t208 * t224) - g(3) * (t209 * pkin(7) + qJ(1) + (pkin(2) * t212 - qJ(3) * t215) * t207) 0, 0, 0, 0, 0, -g(1) * t195 - g(2) * t196 - g(3) * t201, -g(1) * (t199 * t214 - t206 * t223) - g(2) * (t197 * t214 + t208 * t223) - g(3) * (-t209 * t211 - t214 * t220) 0, 0, 0, 0, 0, -g(1) * (t195 * t213 + t200 * t210) - g(2) * (t196 * t213 + t198 * t210) - g(3) * (t201 * t213 + t210 * t222) -g(1) * (-t195 * t210 + t200 * t213) - g(2) * (-t196 * t210 + t198 * t213) - g(3) * (-t201 * t210 + t213 * t222) 0, 0, 0, 0, 0, -g(1) * (t195 * t204 + t200 * t203) - g(2) * (t196 * t204 + t198 * t203) - g(3) * (t201 * t204 + t203 * t222) -g(1) * (-t195 * t203 + t200 * t204) - g(2) * (-t196 * t203 + t198 * t204) - g(3) * (-t201 * t203 + t204 * t222);];
U_reg  = t1;
