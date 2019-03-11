% Calculate minimal parameter regressor of potential energy for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:12
% EndTime: 2019-03-09 03:06:12
% DurationCPUTime: 0.09s
% Computational Cost: add. (139->42), mult. (120->58), div. (0->0), fcn. (123->10), ass. (0->33)
t201 = qJ(3) + pkin(10);
t196 = cos(t201);
t222 = pkin(4) * t196;
t194 = sin(t201);
t221 = g(3) * t194;
t220 = qJ(2) + pkin(6);
t202 = qJ(1) + pkin(9);
t195 = sin(t202);
t204 = sin(qJ(5));
t219 = t195 * t204;
t207 = cos(qJ(5));
t218 = t195 * t207;
t197 = cos(t202);
t217 = t197 * t204;
t216 = t197 * t207;
t205 = sin(qJ(3));
t215 = t205 * pkin(3) + t220;
t208 = cos(qJ(3));
t193 = t208 * pkin(3) + pkin(2);
t203 = -qJ(4) - pkin(7);
t206 = sin(qJ(1));
t214 = t206 * pkin(1) + t195 * t193 + t197 * t203;
t213 = g(1) * t197 + g(2) * t195;
t209 = cos(qJ(1));
t212 = -g(1) * t209 - g(2) * t206;
t211 = t209 * pkin(1) + t197 * t193 - t195 * t203;
t185 = t196 * t219 + t216;
t187 = t196 * t217 - t218;
t210 = g(1) * t187 + g(2) * t185 + t204 * t221;
t188 = t196 * t216 + t219;
t186 = t196 * t218 - t217;
t184 = -g(1) * t188 - g(2) * t186 - t207 * t221;
t1 = [0, t212, g(1) * t206 - g(2) * t209, pkin(1) * t212 - g(3) * t220, 0, 0, 0, 0, 0, -g(3) * t205 - t208 * t213, -g(3) * t208 + t205 * t213, -g(1) * t195 + g(2) * t197, -g(1) * t211 - g(2) * t214 - g(3) * t215, 0, 0, 0, 0, 0, t184, t210, t184, g(3) * t196 - t194 * t213, -t210, -g(1) * (t188 * pkin(5) + t187 * qJ(6) + t197 * t222 + t211) - g(2) * (t186 * pkin(5) + t185 * qJ(6) + t195 * t222 + t214) - g(3) * (-t196 * pkin(8) + t215) + (-g(3) * (pkin(5) * t207 + qJ(6) * t204 + pkin(4)) - t213 * pkin(8)) * t194;];
U_reg  = t1;
