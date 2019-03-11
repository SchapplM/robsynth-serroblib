% Calculate minimal parameter regressor of potential energy for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:24
% EndTime: 2019-03-09 02:45:24
% DurationCPUTime: 0.10s
% Computational Cost: add. (102->41), mult. (108->54), div. (0->0), fcn. (103->8), ass. (0->25)
t168 = sin(qJ(3));
t185 = qJ(4) * t168 + pkin(2);
t171 = cos(qJ(3));
t184 = g(3) * t171;
t183 = qJ(2) + pkin(6);
t166 = qJ(1) + pkin(9);
t160 = sin(t166);
t181 = t160 * t171;
t161 = cos(t166);
t180 = t161 * t171;
t167 = sin(qJ(6));
t179 = t167 * t168;
t170 = cos(qJ(6));
t178 = t168 * t170;
t169 = sin(qJ(1));
t177 = t169 * pkin(1) + pkin(3) * t181 + t185 * t160;
t176 = g(1) * t161 + g(2) * t160;
t172 = cos(qJ(1));
t175 = -g(1) * t172 - g(2) * t169;
t174 = t172 * pkin(1) + pkin(3) * t180 + t160 * pkin(7) + t185 * t161;
t173 = t168 * pkin(3) - t171 * qJ(4) + t183;
t152 = g(1) * t160 - g(2) * t161;
t151 = g(3) * t168 + t176 * t171;
t150 = t176 * t168 - t184;
t1 = [0, t175, g(1) * t169 - g(2) * t172, t175 * pkin(1) - g(3) * t183, 0, 0, 0, 0, 0, -t151, t150, -t151, -t152, -t150, -g(1) * t174 - g(2) * (-t161 * pkin(7) + t177) - g(3) * t173, -t150, t151, t152, -g(1) * (pkin(4) * t180 - t160 * qJ(5) + t174) - g(2) * (pkin(4) * t181 + (-pkin(7) + qJ(5)) * t161 + t177) - g(3) * (t168 * pkin(4) + t173) 0, 0, 0, 0, 0, -g(1) * (-t160 * t167 + t161 * t178) - g(2) * (t160 * t178 + t161 * t167) + t170 * t184, -g(1) * (-t160 * t170 - t161 * t179) - g(2) * (-t160 * t179 + t161 * t170) - t167 * t184;];
U_reg  = t1;
