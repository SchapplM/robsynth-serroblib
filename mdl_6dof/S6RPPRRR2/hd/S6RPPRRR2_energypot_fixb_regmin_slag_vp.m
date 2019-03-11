% Calculate minimal parameter regressor of potential energy for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:21
% EndTime: 2019-03-09 02:21:21
% DurationCPUTime: 0.07s
% Computational Cost: add. (95->35), mult. (80->54), div. (0->0), fcn. (82->12), ass. (0->28)
t193 = pkin(11) + qJ(4);
t187 = sin(t193);
t214 = g(3) * t187;
t213 = g(3) * (qJ(2) + pkin(6));
t194 = qJ(1) + pkin(10);
t188 = sin(t194);
t195 = qJ(5) + qJ(6);
t191 = sin(t195);
t212 = t188 * t191;
t192 = cos(t195);
t211 = t188 * t192;
t199 = sin(qJ(5));
t210 = t188 * t199;
t201 = cos(qJ(5));
t209 = t188 * t201;
t190 = cos(t194);
t208 = t190 * t191;
t207 = t190 * t192;
t206 = t190 * t199;
t205 = t190 * t201;
t204 = g(1) * t190 + g(2) * t188;
t200 = sin(qJ(1));
t202 = cos(qJ(1));
t203 = -g(1) * t202 - g(2) * t200;
t197 = cos(pkin(11));
t196 = sin(pkin(11));
t189 = cos(t193);
t1 = [0, t203, g(1) * t200 - g(2) * t202, t203 * pkin(1) - t213, -g(3) * t196 - t204 * t197, -g(3) * t197 + t204 * t196, -g(1) * t188 + g(2) * t190, -g(1) * (t202 * pkin(1) + t190 * pkin(2) + t188 * qJ(3)) - g(2) * (t200 * pkin(1) + t188 * pkin(2) - t190 * qJ(3)) - t213, 0, 0, 0, 0, 0, -t204 * t189 - t214, -g(3) * t189 + t204 * t187, 0, 0, 0, 0, 0, -g(1) * (t189 * t205 + t210) - g(2) * (t189 * t209 - t206) - t201 * t214, -g(1) * (-t189 * t206 + t209) - g(2) * (-t189 * t210 - t205) + t199 * t214, 0, 0, 0, 0, 0, -g(1) * (t189 * t207 + t212) - g(2) * (t189 * t211 - t208) - t192 * t214, -g(1) * (-t189 * t208 + t211) - g(2) * (-t189 * t212 - t207) + t191 * t214;];
U_reg  = t1;
