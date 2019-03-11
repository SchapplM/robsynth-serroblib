% Calculate minimal parameter regressor of potential energy for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:57:25
% EndTime: 2019-03-09 02:57:25
% DurationCPUTime: 0.08s
% Computational Cost: add. (74->42), mult. (97->53), div. (0->0), fcn. (92->8), ass. (0->26)
t173 = qJ(3) + pkin(9);
t167 = sin(t173);
t168 = cos(t173);
t194 = pkin(4) * t167 - qJ(5) * t168;
t176 = sin(qJ(3));
t193 = pkin(3) * t176;
t191 = g(3) * t167;
t175 = sin(qJ(6));
t177 = sin(qJ(1));
t189 = t177 * t175;
t178 = cos(qJ(6));
t188 = t177 * t178;
t180 = cos(qJ(1));
t187 = t180 * t175;
t186 = t180 * t178;
t185 = t180 * pkin(1) + t177 * qJ(2);
t179 = cos(qJ(3));
t184 = t179 * pkin(3) + pkin(2) + pkin(6);
t183 = t177 * t193 + t185;
t182 = -qJ(2) - t193;
t170 = t177 * pkin(1);
t174 = -qJ(4) - pkin(7);
t181 = -t177 * t174 + t170;
t164 = g(1) * t177 - g(2) * t180;
t165 = g(1) * t180 + g(2) * t177;
t1 = [0, -t165, t164, t165, -t164, -g(1) * t185 - g(2) * (-t180 * qJ(2) + t170) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t179 - t164 * t176, g(3) * t176 - t164 * t179, -t165, -g(1) * (-t180 * t174 + t183) - g(2) * (t182 * t180 + t181) - g(3) * t184, -t165, g(3) * t168 + t164 * t167, t164 * t168 - t191, -g(1) * (t194 * t177 + t183) - g(2) * t181 - g(3) * (t168 * pkin(4) + t167 * qJ(5) + t184) + (g(1) * t174 - g(2) * (t182 - t194)) * t180, 0, 0, 0, 0, 0, -g(1) * (-t168 * t189 + t186) - g(2) * (t168 * t187 + t188) - t175 * t191, -g(1) * (-t168 * t188 - t187) - g(2) * (t168 * t186 - t189) - t178 * t191;];
U_reg  = t1;
