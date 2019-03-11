% Calculate minimal parameter regressor of potential energy for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:04:05
% EndTime: 2019-03-09 18:04:05
% DurationCPUTime: 0.05s
% Computational Cost: add. (82->29), mult. (69->41), div. (0->0), fcn. (70->10), ass. (0->22)
t215 = qJ(2) + qJ(3);
t211 = pkin(11) + qJ(5) + t215;
t209 = sin(t211);
t227 = g(3) * t209;
t216 = sin(qJ(6));
t218 = sin(qJ(1));
t226 = t218 * t216;
t219 = cos(qJ(6));
t225 = t218 * t219;
t221 = cos(qJ(1));
t224 = t221 * t216;
t223 = t221 * t219;
t222 = g(1) * t221 + g(2) * t218;
t220 = cos(qJ(2));
t217 = sin(qJ(2));
t214 = -qJ(4) - pkin(8) - pkin(7);
t213 = cos(t215);
t212 = sin(t215);
t210 = cos(t211);
t208 = g(1) * t218 - g(2) * t221;
t207 = t220 * pkin(2) + pkin(3) * t213 + pkin(1);
t1 = [0, -t222, t208, 0, 0, 0, 0, 0, -g(3) * t217 - t222 * t220, -g(3) * t220 + t222 * t217, 0, 0, 0, 0, 0, -g(3) * t212 - t222 * t213, -g(3) * t213 + t222 * t212, -t208, -g(1) * (t221 * t207 - t218 * t214) - g(2) * (t218 * t207 + t221 * t214) - g(3) * (t217 * pkin(2) + pkin(3) * t212 + pkin(6)) 0, 0, 0, 0, 0, -t222 * t210 - t227, -g(3) * t210 + t222 * t209, 0, 0, 0, 0, 0, -g(1) * (t210 * t223 + t226) - g(2) * (t210 * t225 - t224) - t219 * t227, -g(1) * (-t210 * t224 + t225) - g(2) * (-t210 * t226 - t223) + t216 * t227;];
U_reg  = t1;
