% Calculate minimal parameter regressor of potential energy for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:54:37
% EndTime: 2019-03-09 02:54:37
% DurationCPUTime: 0.10s
% Computational Cost: add. (92->49), mult. (110->66), div. (0->0), fcn. (109->10), ass. (0->33)
t192 = qJ(3) + pkin(9);
t184 = sin(t192);
t186 = cos(t192);
t217 = pkin(4) * t184 - qJ(5) * t186;
t196 = sin(qJ(3));
t216 = pkin(3) * t196;
t214 = g(3) * t186;
t191 = pkin(10) + qJ(6);
t183 = sin(t191);
t197 = sin(qJ(1));
t212 = t197 * t183;
t185 = cos(t191);
t211 = t197 * t185;
t193 = sin(pkin(10));
t210 = t197 * t193;
t194 = cos(pkin(10));
t209 = t197 * t194;
t199 = cos(qJ(1));
t208 = t199 * t183;
t207 = t199 * t185;
t206 = t199 * t193;
t205 = t199 * t194;
t204 = t199 * pkin(1) + t197 * qJ(2);
t198 = cos(qJ(3));
t203 = t198 * pkin(3) + pkin(2) + pkin(6);
t202 = t197 * t216 + t204;
t201 = -qJ(2) - t216;
t188 = t197 * pkin(1);
t195 = -qJ(4) - pkin(7);
t200 = -t197 * t195 + t188;
t180 = g(1) * t197 - g(2) * t199;
t181 = g(1) * t199 + g(2) * t197;
t1 = [0, -t181, t180, t181, -t180, -g(1) * t204 - g(2) * (-t199 * qJ(2) + t188) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t198 - t180 * t196, g(3) * t196 - t180 * t198, -t181, -g(1) * (-t199 * t195 + t202) - g(2) * (t201 * t199 + t200) - g(3) * t203, -g(1) * (t184 * t209 + t206) - g(2) * (-t184 * t205 + t210) - t194 * t214, -g(1) * (-t184 * t210 + t205) - g(2) * (t184 * t206 + t209) + t193 * t214, -g(3) * t184 + t180 * t186, -g(1) * (t217 * t197 + t202) - g(2) * t200 - g(3) * (t186 * pkin(4) + t184 * qJ(5) + t203) + (g(1) * t195 - g(2) * (t201 - t217)) * t199, 0, 0, 0, 0, 0, -g(1) * (t184 * t211 + t208) - g(2) * (-t184 * t207 + t212) - t185 * t214, -g(1) * (-t184 * t212 + t207) - g(2) * (t184 * t208 + t211) + t183 * t214;];
U_reg  = t1;
