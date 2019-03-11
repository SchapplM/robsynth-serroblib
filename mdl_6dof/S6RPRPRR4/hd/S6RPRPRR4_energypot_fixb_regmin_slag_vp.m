% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:26
% EndTime: 2019-03-09 03:46:26
% DurationCPUTime: 0.07s
% Computational Cost: add. (87->37), mult. (90->57), div. (0->0), fcn. (92->10), ass. (0->25)
t195 = cos(qJ(3));
t206 = g(3) * t195;
t205 = qJ(2) + pkin(6);
t190 = qJ(5) + qJ(6);
t187 = sin(t190);
t192 = sin(qJ(3));
t204 = t187 * t192;
t188 = cos(t190);
t203 = t188 * t192;
t191 = sin(qJ(5));
t202 = t191 * t192;
t194 = cos(qJ(5));
t201 = t192 * t194;
t189 = qJ(1) + pkin(10);
t185 = sin(t189);
t186 = cos(t189);
t200 = g(1) * t186 + g(2) * t185;
t193 = sin(qJ(1));
t196 = cos(qJ(1));
t199 = -g(1) * t196 - g(2) * t193;
t198 = pkin(3) * t195 + qJ(4) * t192 + pkin(2);
t197 = t199 * pkin(1);
t184 = g(3) * t192 + t200 * t195;
t183 = t200 * t192 - t206;
t1 = [0, t199, g(1) * t193 - g(2) * t196, -g(3) * t205 + t197, 0, 0, 0, 0, 0, -t184, t183, -g(1) * t185 + g(2) * t186, t184, -t183, -g(3) * (t192 * pkin(3) - t195 * qJ(4) + t205) + t197 + (g(2) * pkin(7) - g(1) * t198) * t186 + (-g(1) * pkin(7) - g(2) * t198) * t185, 0, 0, 0, 0, 0, -g(1) * (t185 * t194 + t186 * t202) - g(2) * (t185 * t202 - t186 * t194) + t191 * t206, -g(1) * (-t185 * t191 + t186 * t201) - g(2) * (t185 * t201 + t186 * t191) + t194 * t206, 0, 0, 0, 0, 0, -g(1) * (t185 * t188 + t186 * t204) - g(2) * (t185 * t204 - t186 * t188) + t187 * t206, -g(1) * (-t185 * t187 + t186 * t203) - g(2) * (t185 * t203 + t186 * t187) + t188 * t206;];
U_reg  = t1;
