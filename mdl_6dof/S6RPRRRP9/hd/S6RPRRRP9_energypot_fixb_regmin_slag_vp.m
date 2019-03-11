% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP9
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
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:30
% EndTime: 2019-03-09 06:28:30
% DurationCPUTime: 0.07s
% Computational Cost: add. (74->44), mult. (101->59), div. (0->0), fcn. (103->8), ass. (0->28)
t174 = cos(qJ(3));
t187 = g(3) * t174;
t169 = qJ(4) + qJ(5);
t163 = sin(t169);
t170 = sin(qJ(4));
t186 = t170 * pkin(4) + pkin(5) * t163 + pkin(7);
t172 = sin(qJ(1));
t185 = t172 * t163;
t164 = cos(t169);
t184 = t172 * t164;
t183 = t172 * t170;
t173 = cos(qJ(4));
t182 = t172 * t173;
t175 = cos(qJ(1));
t181 = t175 * t163;
t180 = t175 * t164;
t179 = t175 * t170;
t178 = t175 * t173;
t177 = g(1) * (t175 * pkin(1) + t172 * qJ(2));
t161 = g(1) * t172 - g(2) * t175;
t159 = t173 * pkin(4) + pkin(5) * t164 + pkin(3);
t168 = -qJ(6) - pkin(9) - pkin(8);
t171 = sin(qJ(3));
t176 = t159 * t171 + t168 * t174;
t166 = t172 * pkin(1);
t162 = g(1) * t175 + g(2) * t172;
t158 = -g(3) * t171 + t161 * t174;
t1 = [0, -t162, t161, t162, -t161, -t177 - g(2) * (-t175 * qJ(2) + t166) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t161 * t171 - t187, -t158, 0, 0, 0, 0, 0, -g(1) * (t171 * t182 + t179) - g(2) * (-t171 * t178 + t183) - t173 * t187, -g(1) * (-t171 * t183 + t178) - g(2) * (t171 * t179 + t182) + t170 * t187, 0, 0, 0, 0, 0, -g(1) * (t171 * t184 + t181) - g(2) * (-t171 * t180 + t185) - t164 * t187, -g(1) * (-t171 * t185 + t180) - g(2) * (t171 * t181 + t184) + t163 * t187, t158, -t177 - g(2) * t166 - g(3) * (t174 * t159 - t171 * t168 + pkin(2) + pkin(6)) + (-g(1) * t176 - g(2) * t186) * t172 + (-g(1) * t186 - g(2) * (-qJ(2) - t176)) * t175;];
U_reg  = t1;
