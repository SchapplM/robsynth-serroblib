% Calculate minimal parameter regressor of potential energy for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x38]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:29:47
% EndTime: 2019-03-10 03:29:47
% DurationCPUTime: 0.05s
% Computational Cost: add. (80->22), mult. (64->34), div. (0->0), fcn. (68->12), ass. (0->22)
t211 = qJ(2) + qJ(3);
t210 = qJ(4) + t211;
t207 = qJ(5) + t210;
t203 = sin(t207);
t223 = g(3) * t203;
t212 = sin(qJ(6));
t214 = sin(qJ(1));
t222 = t214 * t212;
t215 = cos(qJ(6));
t221 = t214 * t215;
t217 = cos(qJ(1));
t220 = t217 * t212;
t219 = t217 * t215;
t218 = g(1) * t217 + g(2) * t214;
t216 = cos(qJ(2));
t213 = sin(qJ(2));
t209 = cos(t211);
t208 = sin(t211);
t206 = cos(t210);
t205 = sin(t210);
t204 = cos(t207);
t1 = [0, -t218, g(1) * t214 - g(2) * t217, 0, 0, 0, 0, 0, -g(3) * t213 - t218 * t216, -g(3) * t216 + t218 * t213, 0, 0, 0, 0, 0, -g(3) * t208 - t218 * t209, -g(3) * t209 + t218 * t208, 0, 0, 0, 0, 0, -g(3) * t205 - t218 * t206, -g(3) * t206 + t218 * t205, 0, 0, 0, 0, 0, -t218 * t204 - t223, -g(3) * t204 + t218 * t203, 0, 0, 0, 0, 0, -g(1) * (t204 * t219 + t222) - g(2) * (t204 * t221 - t220) - t215 * t223, -g(1) * (-t204 * t220 + t221) - g(2) * (-t204 * t222 - t219) + t212 * t223;];
U_reg  = t1;
