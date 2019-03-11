% Calculate minimal parameter regressor of potential energy for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:47
% EndTime: 2019-03-08 19:53:47
% DurationCPUTime: 0.12s
% Computational Cost: add. (127->61), mult. (296->92), div. (0->0), fcn. (359->10), ass. (0->35)
t202 = sin(pkin(6));
t224 = pkin(7) * t202;
t206 = sin(qJ(4));
t223 = t202 * t206;
t207 = sin(qJ(2));
t222 = t202 * t207;
t209 = cos(qJ(4));
t221 = t202 * t209;
t210 = cos(qJ(2));
t220 = t202 * t210;
t204 = cos(pkin(6));
t219 = t204 * t207;
t218 = t204 * t210;
t217 = pkin(2) * t222 + t204 * pkin(7) + qJ(1);
t201 = sin(pkin(10));
t203 = cos(pkin(10));
t188 = t201 * t207 - t203 * t218;
t189 = t201 * t210 + t203 * t219;
t216 = t201 * pkin(1) + t189 * pkin(2) + t188 * qJ(3);
t190 = t201 * t218 + t203 * t207;
t191 = -t201 * t219 + t203 * t210;
t215 = t203 * pkin(1) + t191 * pkin(2) + t190 * qJ(3) + t201 * t224;
t180 = -t190 * t209 + t201 * t223;
t182 = t188 * t209 + t203 * t223;
t192 = t204 * t206 + t209 * t220;
t214 = g(1) * t180 - g(2) * t182 + g(3) * t192;
t181 = t190 * t206 + t201 * t221;
t183 = -t188 * t206 + t203 * t221;
t193 = t204 * t209 - t206 * t220;
t213 = g(1) * t181 - g(2) * t183 + g(3) * t193;
t212 = -g(1) * t190 - g(2) * t188 + g(3) * t220;
t211 = g(1) * t191 + g(2) * t189 + g(3) * t222;
t208 = cos(qJ(6));
t205 = sin(qJ(6));
t1 = [-g(3) * qJ(1), 0, -t211, -t212, t211, t212, -g(1) * t215 - g(2) * (-t203 * t224 + t216) - g(3) * (-qJ(3) * t220 + t217) 0, 0, 0, 0, 0, -t213, t214, -t211, t213, -t214, -g(1) * (t181 * pkin(4) + t191 * pkin(8) + t180 * qJ(5) + t215) - g(2) * (-t183 * pkin(4) + t189 * pkin(8) - t182 * qJ(5) + t216) - g(3) * (t204 * pkin(3) + t193 * pkin(4) + t192 * qJ(5) + t217) + (-g(1) * pkin(3) * t201 - g(3) * (pkin(8) * t207 - qJ(3) * t210) - g(2) * (-pkin(3) - pkin(7)) * t203) * t202, 0, 0, 0, 0, 0, -g(1) * (t180 * t205 + t191 * t208) - g(2) * (-t182 * t205 + t189 * t208) - g(3) * (t192 * t205 + t208 * t222) -g(1) * (t180 * t208 - t191 * t205) - g(2) * (-t182 * t208 - t189 * t205) - g(3) * (t192 * t208 - t205 * t222);];
U_reg  = t1;
