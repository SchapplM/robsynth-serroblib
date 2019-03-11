% Calculate minimal parameter regressor of potential energy for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:35
% EndTime: 2019-03-09 12:48:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (85->52), mult. (127->68), div. (0->0), fcn. (129->8), ass. (0->31)
t212 = cos(qJ(2));
t228 = g(3) * t212;
t209 = sin(qJ(2));
t227 = t209 * pkin(2) + pkin(6);
t210 = sin(qJ(1));
t226 = t209 * t210;
t207 = qJ(4) + qJ(5);
t200 = sin(t207);
t225 = t210 * t200;
t201 = cos(t207);
t224 = t210 * t201;
t208 = sin(qJ(4));
t223 = t210 * t208;
t211 = cos(qJ(4));
t222 = t210 * t211;
t221 = t210 * t212;
t213 = cos(qJ(1));
t220 = t213 * t200;
t219 = t213 * t201;
t218 = t213 * t208;
t217 = t213 * t211;
t216 = t210 * pkin(1) + pkin(2) * t221 + qJ(3) * t226;
t215 = t210 * pkin(7) + (pkin(2) * t212 + qJ(3) * t209 + pkin(1)) * t213;
t214 = g(1) * t213 + g(2) * t210;
t206 = -qJ(6) - pkin(9) - pkin(8);
t195 = g(1) * t210 - g(2) * t213;
t194 = t208 * pkin(4) + pkin(5) * t200;
t193 = t211 * pkin(4) + pkin(5) * t201 + pkin(3);
t192 = g(3) * t209 + t214 * t212;
t191 = t214 * t209 - t228;
t1 = [0, -t214, t195, 0, 0, 0, 0, 0, -t192, t191, -t195, t192, -t191, -g(1) * t215 - g(2) * (-t213 * pkin(7) + t216) - g(3) * (-t212 * qJ(3) + t227) 0, 0, 0, 0, 0, -g(1) * (t209 * t218 + t222) - g(2) * (t209 * t223 - t217) + t208 * t228, -g(1) * (t209 * t217 - t223) - g(2) * (t209 * t222 + t218) + t211 * t228, 0, 0, 0, 0, 0, -g(1) * (t209 * t220 + t224) - g(2) * (t209 * t225 - t219) + t200 * t228, -g(1) * (t209 * t219 - t225) - g(2) * (t209 * t224 + t220) + t201 * t228, -t192, -g(1) * (t210 * t193 + t215) - g(2) * (t194 * t226 - t206 * t221 + t216) - g(3) * (-t209 * t206 + (-qJ(3) - t194) * t212 + t227) + (-g(1) * (t194 * t209 - t206 * t212) - g(2) * (-pkin(7) - t193)) * t213;];
U_reg  = t1;
