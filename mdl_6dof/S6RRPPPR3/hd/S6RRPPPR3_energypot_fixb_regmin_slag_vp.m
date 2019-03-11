% Calculate minimal parameter regressor of potential energy for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:38
% EndTime: 2019-03-09 08:15:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (90->53), mult. (158->67), div. (0->0), fcn. (157->8), ass. (0->34)
t200 = sin(qJ(2));
t224 = qJ(3) * t200 + pkin(1);
t202 = cos(qJ(2));
t223 = g(3) * t202;
t222 = t200 * pkin(2) + pkin(6);
t197 = pkin(9) + qJ(6);
t188 = sin(t197);
t201 = sin(qJ(1));
t220 = t201 * t188;
t189 = cos(t197);
t219 = t201 * t189;
t198 = sin(pkin(9));
t218 = t201 * t198;
t199 = cos(pkin(9));
t217 = t201 * t199;
t216 = t201 * t202;
t203 = cos(qJ(1));
t215 = t202 * t203;
t214 = t203 * t188;
t213 = t203 * t189;
t212 = t203 * t198;
t211 = t203 * t199;
t210 = pkin(2) * t215 + t201 * pkin(7) + t224 * t203;
t209 = -t202 * qJ(3) + t222;
t208 = g(1) * t203 + g(2) * t201;
t207 = pkin(4) * t200 + qJ(5) * t202;
t206 = pkin(2) * t216 - t203 * pkin(7) + t224 * t201;
t205 = pkin(3) * t216 + t203 * qJ(4) + t206;
t204 = pkin(3) * t215 - t201 * qJ(4) + t210;
t191 = t200 * pkin(3);
t181 = g(1) * t201 - g(2) * t203;
t180 = g(3) * t200 + t208 * t202;
t179 = t208 * t200 - t223;
t1 = [0, -t208, t181, 0, 0, 0, 0, 0, -t180, t179, -t180, -t181, -t179, -g(1) * t210 - g(2) * t206 - g(3) * t209, -t179, t180, t181, -g(1) * t204 - g(2) * t205 - g(3) * (t191 + t209) -g(1) * (t200 * t211 - t218) - g(2) * (t200 * t217 + t212) + t199 * t223, -g(1) * (-t200 * t212 - t217) - g(2) * (-t200 * t218 + t211) - t198 * t223, -t180, -g(1) * (t207 * t203 + t204) - g(2) * (t207 * t201 + t205) - g(3) * (t200 * qJ(5) + t191 + (-pkin(4) - qJ(3)) * t202 + t222) 0, 0, 0, 0, 0, -g(1) * (t200 * t213 - t220) - g(2) * (t200 * t219 + t214) + t189 * t223, -g(1) * (-t200 * t214 - t219) - g(2) * (-t200 * t220 + t213) - t188 * t223;];
U_reg  = t1;
