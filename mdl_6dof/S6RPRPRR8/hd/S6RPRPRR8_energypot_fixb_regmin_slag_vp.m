% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:30
% EndTime: 2019-03-09 03:59:30
% DurationCPUTime: 0.08s
% Computational Cost: add. (63->39), mult. (81->53), div. (0->0), fcn. (83->10), ass. (0->27)
t195 = sin(qJ(3));
t210 = pkin(3) * t195;
t191 = qJ(3) + pkin(10);
t209 = g(3) * cos(t191);
t192 = qJ(5) + qJ(6);
t186 = sin(t192);
t199 = cos(qJ(1));
t208 = t186 * t199;
t187 = cos(t192);
t207 = t187 * t199;
t194 = sin(qJ(5));
t206 = t194 * t199;
t196 = sin(qJ(1));
t205 = t196 * t186;
t204 = t196 * t187;
t203 = t196 * t194;
t197 = cos(qJ(5));
t202 = t196 * t197;
t201 = t197 * t199;
t200 = t199 * pkin(1) + t196 * qJ(2);
t182 = g(1) * t196 - g(2) * t199;
t198 = cos(qJ(3));
t193 = -qJ(4) - pkin(7);
t189 = t196 * pkin(1);
t184 = sin(t191);
t183 = g(1) * t199 + g(2) * t196;
t1 = [0, -t183, t182, t183, -t182, -g(1) * t200 - g(2) * (-qJ(2) * t199 + t189) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t198 - t182 * t195, g(3) * t195 - t182 * t198, -t183, -g(1) * (-t193 * t199 + t196 * t210 + t200) - g(2) * (-t196 * t193 + t189 + (-qJ(2) - t210) * t199) - g(3) * (pkin(3) * t198 + pkin(2) + pkin(6)) 0, 0, 0, 0, 0, -g(1) * (t184 * t202 + t206) - g(2) * (-t184 * t201 + t203) - t197 * t209, -g(1) * (-t184 * t203 + t201) - g(2) * (t184 * t206 + t202) + t194 * t209, 0, 0, 0, 0, 0, -g(1) * (t184 * t204 + t208) - g(2) * (-t184 * t207 + t205) - t187 * t209, -g(1) * (-t184 * t205 + t207) - g(2) * (t184 * t208 + t204) + t186 * t209;];
U_reg  = t1;
