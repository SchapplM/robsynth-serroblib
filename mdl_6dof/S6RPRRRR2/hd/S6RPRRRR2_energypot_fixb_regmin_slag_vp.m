% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:58:29
% EndTime: 2019-03-09 06:58:30
% DurationCPUTime: 0.06s
% Computational Cost: add. (81->28), mult. (69->48), div. (0->0), fcn. (74->12), ass. (0->23)
t197 = qJ(3) + qJ(4);
t192 = sin(t197);
t210 = g(3) * t192;
t196 = qJ(5) + qJ(6);
t191 = sin(t196);
t194 = cos(t197);
t209 = t191 * t194;
t193 = cos(t196);
t208 = t193 * t194;
t198 = sin(qJ(5));
t207 = t194 * t198;
t201 = cos(qJ(5));
t206 = t194 * t201;
t195 = qJ(1) + pkin(11);
t189 = sin(t195);
t190 = cos(t195);
t205 = g(1) * t190 + g(2) * t189;
t200 = sin(qJ(1));
t203 = cos(qJ(1));
t204 = -g(1) * t203 - g(2) * t200;
t202 = cos(qJ(3));
t199 = sin(qJ(3));
t1 = [0, t204, g(1) * t200 - g(2) * t203, -g(3) * (qJ(2) + pkin(6)) + t204 * pkin(1), 0, 0, 0, 0, 0, -g(3) * t199 - t205 * t202, -g(3) * t202 + t205 * t199, 0, 0, 0, 0, 0, -t205 * t194 - t210, -g(3) * t194 + t205 * t192, 0, 0, 0, 0, 0, -g(1) * (t189 * t198 + t190 * t206) - g(2) * (t189 * t206 - t190 * t198) - t201 * t210, -g(1) * (t189 * t201 - t190 * t207) - g(2) * (-t189 * t207 - t190 * t201) + t198 * t210, 0, 0, 0, 0, 0, -g(1) * (t189 * t191 + t190 * t208) - g(2) * (t189 * t208 - t190 * t191) - t193 * t210, -g(1) * (t189 * t193 - t190 * t209) - g(2) * (-t189 * t209 - t190 * t193) + t191 * t210;];
U_reg  = t1;
