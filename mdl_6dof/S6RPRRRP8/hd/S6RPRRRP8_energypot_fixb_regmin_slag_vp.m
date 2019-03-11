% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:33
% EndTime: 2019-03-09 06:24:33
% DurationCPUTime: 0.11s
% Computational Cost: add. (95->44), mult. (124->55), div. (0->0), fcn. (130->8), ass. (0->29)
t173 = qJ(3) + qJ(4);
t168 = sin(t173);
t175 = sin(qJ(3));
t191 = pkin(3) * t175 + pkin(4) * t168;
t176 = sin(qJ(1));
t188 = g(1) * t176;
t169 = cos(t173);
t187 = g(3) * t169;
t174 = sin(qJ(5));
t186 = t176 * t174;
t177 = cos(qJ(5));
t185 = t176 * t177;
t179 = cos(qJ(1));
t184 = t179 * t174;
t183 = t179 * t177;
t182 = t179 * pkin(1) + t176 * qJ(2);
t165 = -g(2) * t179 + t188;
t161 = t168 * t186 - t183;
t163 = t168 * t184 + t185;
t181 = g(1) * t161 - g(2) * t163 + t174 * t187;
t180 = -pkin(8) - pkin(7);
t178 = cos(qJ(3));
t171 = t176 * pkin(1);
t166 = g(1) * t179 + g(2) * t176;
t164 = -t168 * t183 + t186;
t162 = t168 * t185 + t184;
t160 = -g(3) * t168 + t165 * t169;
t159 = -g(1) * t162 - g(2) * t164 - t177 * t187;
t1 = [0, -t166, t165, t166, -t165, -g(1) * t182 - g(2) * (-t179 * qJ(2) + t171) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t178 - t165 * t175, g(3) * t175 - t165 * t178, 0, 0, 0, 0, 0, -t165 * t168 - t187, -t160, 0, 0, 0, 0, 0, t159, t181, t159, t160, -t181, -g(1) * (t162 * pkin(5) + t161 * qJ(6) + t191 * t176 + t182) - g(2) * (t164 * pkin(5) - t163 * qJ(6) - t176 * t180 + t171) - g(3) * (t178 * pkin(3) + t168 * pkin(9) + pkin(2) + pkin(6)) + (pkin(9) * t188 - g(3) * (pkin(5) * t177 + qJ(6) * t174 + pkin(4))) * t169 + (g(1) * t180 - g(2) * (pkin(9) * t169 - qJ(2) - t191)) * t179;];
U_reg  = t1;
