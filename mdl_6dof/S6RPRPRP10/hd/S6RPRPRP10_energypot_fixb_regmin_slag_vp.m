% Calculate minimal parameter regressor of potential energy for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:31
% EndTime: 2019-03-09 03:32:31
% DurationCPUTime: 0.08s
% Computational Cost: add. (76->51), mult. (147->61), div. (0->0), fcn. (150->6), ass. (0->30)
t169 = sin(qJ(1));
t171 = cos(qJ(3));
t172 = cos(qJ(1));
t185 = t172 * t171 * qJ(4) + t169 * pkin(7);
t184 = g(1) * t169;
t183 = g(2) * t172;
t168 = sin(qJ(3));
t182 = g(3) * t168;
t181 = t168 * t169;
t180 = t169 * t171;
t167 = sin(qJ(5));
t179 = t172 * t167;
t170 = cos(qJ(5));
t178 = t172 * t170;
t177 = t172 * pkin(1) + t169 * qJ(2);
t176 = t171 * pkin(3) + t168 * qJ(4) + pkin(2) + pkin(6);
t163 = t169 * pkin(1);
t175 = -t172 * qJ(2) + t163;
t154 = -t183 + t184;
t174 = pkin(3) * t181 + t172 * pkin(7) - qJ(4) * t180 + t177;
t150 = t169 * t167 - t171 * t178;
t152 = t170 * t180 + t179;
t173 = -g(1) * t152 - g(2) * t150 + t170 * t182;
t155 = g(1) * t172 + g(2) * t169;
t153 = -t167 * t180 + t178;
t151 = t169 * t170 + t171 * t179;
t149 = t154 * t171 - t182;
t148 = g(1) * t181 + g(3) * t171 - t168 * t183;
t147 = -g(1) * t153 - g(2) * t151 - t167 * t182;
t1 = [0, -t155, t154, t155, -t154, -g(3) * pkin(6) - g(1) * t177 - g(2) * t175, 0, 0, 0, 0, 0, -t148, -t149, -t155, t148, t149, -g(1) * t174 - g(2) * (t163 + (-pkin(3) * t168 - qJ(2)) * t172 + t185) - g(3) * t176, 0, 0, 0, 0, 0, t147, -t173, t147, -t148, t173, -g(1) * (t172 * pkin(4) + t153 * pkin(5) + t152 * qJ(6) + t174) - g(2) * (t169 * pkin(4) + t151 * pkin(5) + t150 * qJ(6) + t175 + t185) - g(3) * (t171 * pkin(8) + t176) + (-pkin(8) * t184 - g(3) * (pkin(5) * t167 - qJ(6) * t170) - (-pkin(3) - pkin(8)) * t183) * t168;];
U_reg  = t1;
