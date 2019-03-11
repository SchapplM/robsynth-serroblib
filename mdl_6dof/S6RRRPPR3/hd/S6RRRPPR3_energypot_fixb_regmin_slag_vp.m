% Calculate minimal parameter regressor of potential energy for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:30:07
% EndTime: 2019-03-09 15:30:07
% DurationCPUTime: 0.08s
% Computational Cost: add. (102->41), mult. (115->50), div. (0->0), fcn. (113->8), ass. (0->26)
t194 = qJ(2) + qJ(3);
t191 = sin(t194);
t199 = cos(qJ(2));
t214 = t199 * pkin(2) + qJ(4) * t191 + pkin(1);
t192 = cos(t194);
t213 = g(3) * t192;
t197 = sin(qJ(1));
t211 = t192 * t197;
t200 = cos(qJ(1));
t210 = t192 * t200;
t195 = sin(qJ(6));
t209 = t197 * t195;
t198 = cos(qJ(6));
t208 = t197 * t198;
t207 = t200 * t195;
t206 = t200 * t198;
t205 = pkin(3) * t210 + t214 * t200;
t201 = -pkin(8) - pkin(7);
t204 = pkin(3) * t211 + t214 * t197 + t200 * t201;
t203 = g(1) * t200 + g(2) * t197;
t196 = sin(qJ(2));
t202 = t196 * pkin(2) + t191 * pkin(3) - t192 * qJ(4) + pkin(6);
t180 = g(1) * t197 - g(2) * t200;
t179 = g(3) * t191 + t203 * t192;
t178 = t203 * t191 - t213;
t1 = [0, -t203, t180, 0, 0, 0, 0, 0, -g(3) * t196 - t203 * t199, -g(3) * t199 + t203 * t196, 0, 0, 0, 0, 0, -t179, t178, -t179, -t180, -t178, -g(1) * (-t197 * t201 + t205) - g(2) * t204 - g(3) * t202, -t178, t179, t180, -g(1) * (pkin(4) * t210 + (-qJ(5) - t201) * t197 + t205) - g(2) * (pkin(4) * t211 + t200 * qJ(5) + t204) - g(3) * (t191 * pkin(4) + t202) 0, 0, 0, 0, 0, -g(1) * (t191 * t206 - t209) - g(2) * (t191 * t208 + t207) + t198 * t213, -g(1) * (-t191 * t207 - t208) - g(2) * (-t191 * t209 + t206) - t195 * t213;];
U_reg  = t1;
