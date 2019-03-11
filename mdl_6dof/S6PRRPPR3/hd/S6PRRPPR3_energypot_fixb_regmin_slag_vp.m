% Calculate minimal parameter regressor of potential energy for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:10:56
% EndTime: 2019-03-08 21:10:56
% DurationCPUTime: 0.12s
% Computational Cost: add. (158->61), mult. (371->87), div. (0->0), fcn. (457->10), ass. (0->35)
t203 = sin(pkin(6));
t225 = pkin(7) * t203;
t224 = pkin(8) - qJ(5);
t223 = cos(pkin(6));
t206 = sin(qJ(3));
t222 = t203 * t206;
t207 = sin(qJ(2));
t221 = t203 * t207;
t209 = cos(qJ(3));
t220 = t203 * t209;
t210 = cos(qJ(2));
t219 = t203 * t210;
t218 = t207 * t223;
t217 = t210 * t223;
t192 = t206 * t221 - t223 * t209;
t193 = t223 * t206 + t207 * t220;
t216 = pkin(2) * t221 + t193 * pkin(3) + t223 * pkin(7) + t192 * qJ(4) + qJ(1);
t202 = sin(pkin(10));
t204 = cos(pkin(10));
t191 = -t202 * t218 + t204 * t210;
t183 = t191 * t206 - t202 * t220;
t184 = t191 * t209 + t202 * t222;
t215 = t204 * pkin(1) + t191 * pkin(2) + t184 * pkin(3) + t183 * qJ(4) + t202 * t225;
t189 = t202 * t210 + t204 * t218;
t181 = t189 * t206 + t204 * t220;
t214 = g(1) * t183 + g(2) * t181 + g(3) * t192;
t182 = t189 * t209 - t204 * t222;
t213 = g(1) * t184 + g(2) * t182 + g(3) * t193;
t188 = t202 * t207 - t204 * t217;
t190 = t202 * t217 + t204 * t207;
t212 = -g(1) * t190 - g(2) * t188 + g(3) * t219;
t211 = t202 * pkin(1) + t189 * pkin(2) + t182 * pkin(3) + t181 * qJ(4) - t204 * t225;
t208 = cos(qJ(6));
t205 = sin(qJ(6));
t1 = [-g(3) * qJ(1), 0, -g(1) * t191 - g(2) * t189 - g(3) * t221, -t212, 0, 0, 0, 0, 0, -t213, t214, -t213, t212, -t214, -g(1) * (t190 * pkin(8) + t215) - g(2) * (t188 * pkin(8) + t211) - g(3) * (-pkin(8) * t219 + t216) -t214, t213, -t212, -g(1) * (t184 * pkin(4) + t224 * t190 + t215) - g(2) * (t182 * pkin(4) + t224 * t188 + t211) - g(3) * (t193 * pkin(4) - t224 * t219 + t216) 0, 0, 0, 0, 0, -g(1) * (t183 * t208 - t190 * t205) - g(2) * (t181 * t208 - t188 * t205) - g(3) * (t192 * t208 + t205 * t219) -g(1) * (-t183 * t205 - t190 * t208) - g(2) * (-t181 * t205 - t188 * t208) - g(3) * (-t192 * t205 + t208 * t219);];
U_reg  = t1;
