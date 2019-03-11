% Calculate minimal parameter regressor of potential energy for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:58:13
% EndTime: 2019-03-10 00:58:13
% DurationCPUTime: 0.07s
% Computational Cost: add. (103->37), mult. (89->47), div. (0->0), fcn. (90->10), ass. (0->24)
t200 = qJ(2) + qJ(3);
t198 = qJ(4) + t200;
t193 = sin(t198);
t215 = g(3) * t193;
t202 = sin(qJ(5));
t204 = sin(qJ(1));
t214 = t204 * t202;
t205 = cos(qJ(5));
t213 = t204 * t205;
t207 = cos(qJ(1));
t212 = t207 * t202;
t211 = t207 * t205;
t210 = pkin(5) * t202 + pkin(7) + pkin(8) + pkin(9);
t209 = g(1) * t207 + g(2) * t204;
t194 = cos(t198);
t195 = t205 * pkin(5) + pkin(4);
t197 = cos(t200);
t201 = -qJ(6) - pkin(10);
t206 = cos(qJ(2));
t208 = -t206 * pkin(2) - pkin(3) * t197 + t193 * t201 - t194 * t195 - pkin(1);
t203 = sin(qJ(2));
t196 = sin(t200);
t191 = -g(3) * t194 + t209 * t193;
t1 = [0, -t209, g(1) * t204 - g(2) * t207, 0, 0, 0, 0, 0, -g(3) * t203 - t209 * t206, -g(3) * t206 + t209 * t203, 0, 0, 0, 0, 0, -g(3) * t196 - t209 * t197, -g(3) * t197 + t209 * t196, 0, 0, 0, 0, 0, -t209 * t194 - t215, t191, 0, 0, 0, 0, 0, -g(1) * (t194 * t211 + t214) - g(2) * (t194 * t213 - t212) - t205 * t215, -g(1) * (-t194 * t212 + t213) - g(2) * (-t194 * t214 - t211) + t202 * t215, -t191, -g(3) * (t203 * pkin(2) + pkin(3) * t196 + t193 * t195 + t194 * t201 + pkin(6)) + (g(1) * t208 + g(2) * t210) * t207 + (-g(1) * t210 + g(2) * t208) * t204;];
U_reg  = t1;
