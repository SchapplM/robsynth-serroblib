% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:03:59
% EndTime: 2019-03-09 22:03:59
% DurationCPUTime: 0.06s
% Computational Cost: add. (104->35), mult. (89->45), div. (0->0), fcn. (90->10), ass. (0->24)
t203 = qJ(2) + qJ(3);
t201 = qJ(4) + t203;
t198 = cos(t201);
t216 = g(3) * t198;
t204 = sin(qJ(6));
t206 = sin(qJ(1));
t215 = t206 * t204;
t207 = cos(qJ(6));
t214 = t206 * t207;
t209 = cos(qJ(1));
t213 = t209 * t204;
t212 = t209 * t207;
t211 = g(1) * t209 + g(2) * t206;
t197 = sin(t201);
t200 = cos(t203);
t208 = cos(qJ(2));
t210 = t208 * pkin(2) + pkin(3) * t200 + pkin(4) * t198 + qJ(5) * t197 + pkin(1);
t205 = sin(qJ(2));
t202 = -pkin(9) - pkin(8) - pkin(7);
t199 = sin(t203);
t196 = g(1) * t206 - g(2) * t209;
t194 = g(3) * t197 + t211 * t198;
t193 = t211 * t197 - t216;
t1 = [0, -t211, t196, 0, 0, 0, 0, 0, -g(3) * t205 - t211 * t208, -g(3) * t208 + t211 * t205, 0, 0, 0, 0, 0, -g(3) * t199 - t211 * t200, -g(3) * t200 + t211 * t199, 0, 0, 0, 0, 0, -t194, t193, -t196, t194, -t193, -g(3) * (t205 * pkin(2) + pkin(3) * t199 + t197 * pkin(4) - t198 * qJ(5) + pkin(6)) + (-g(1) * t210 - g(2) * t202) * t209 + (g(1) * t202 - g(2) * t210) * t206, 0, 0, 0, 0, 0, -g(1) * (t197 * t213 + t214) - g(2) * (t197 * t215 - t212) + t204 * t216, -g(1) * (t197 * t212 - t215) - g(2) * (t197 * t214 + t213) + t207 * t216;];
U_reg  = t1;
