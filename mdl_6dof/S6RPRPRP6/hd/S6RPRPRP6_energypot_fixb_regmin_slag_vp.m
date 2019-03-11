% Calculate minimal parameter regressor of potential energy for
% S6RPRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:45
% EndTime: 2019-03-09 03:19:45
% DurationCPUTime: 0.10s
% Computational Cost: add. (112->51), mult. (127->63), div. (0->0), fcn. (122->8), ass. (0->29)
t191 = pkin(9) + qJ(3);
t187 = sin(t191);
t193 = cos(pkin(9));
t213 = t193 * pkin(2) + qJ(4) * t187 + pkin(1);
t196 = sin(qJ(5));
t212 = pkin(5) * t196;
t188 = cos(t191);
t211 = g(3) * t188;
t199 = cos(qJ(1));
t209 = t188 * t199;
t197 = sin(qJ(1));
t208 = t197 * t196;
t198 = cos(qJ(5));
t207 = t197 * t198;
t206 = t199 * t196;
t205 = t199 * t198;
t192 = sin(pkin(9));
t204 = t192 * pkin(2) + t187 * pkin(3) + pkin(6);
t203 = t187 * t206;
t202 = pkin(3) * t209 + t213 * t199;
t195 = -pkin(7) - qJ(2);
t201 = t199 * t195 + (pkin(3) * t188 + t213) * t197;
t200 = g(1) * t199 + g(2) * t197;
t194 = -qJ(6) - pkin(8);
t186 = t198 * pkin(5) + pkin(4);
t180 = g(1) * t197 - g(2) * t199;
t175 = g(3) * t187 + t200 * t188;
t174 = t200 * t187 - t211;
t1 = [0, -t200, t180, -g(3) * t192 - t200 * t193, -g(3) * t193 + t200 * t192, -t180, -g(1) * (t199 * pkin(1) + t197 * qJ(2)) - g(2) * (t197 * pkin(1) - t199 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t175, t174, -t180, t175, -t174, -g(1) * (-t197 * t195 + t202) - g(2) * t201 - g(3) * (-t188 * qJ(4) + t204) 0, 0, 0, 0, 0, -g(1) * (t203 + t207) - g(2) * (t187 * t208 - t205) + t196 * t211, -g(1) * (t187 * t205 - t208) - g(2) * (t187 * t207 + t206) + t198 * t211, -t175, -g(1) * (pkin(5) * t203 - t194 * t209 + t202) - g(2) * (-t199 * t186 + t201) - g(3) * (-t187 * t194 + (-qJ(4) - t212) * t188 + t204) + (-g(1) * (t186 - t195) - g(2) * (t187 * t212 - t188 * t194)) * t197;];
U_reg  = t1;
