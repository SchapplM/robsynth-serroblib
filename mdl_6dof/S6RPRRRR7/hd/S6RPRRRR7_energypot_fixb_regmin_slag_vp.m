% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x34]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:17:54
% EndTime: 2019-03-09 07:17:54
% DurationCPUTime: 0.05s
% Computational Cost: add. (58->24), mult. (65->37), div. (0->0), fcn. (66->10), ass. (0->20)
t174 = qJ(3) + qJ(4);
t173 = qJ(5) + t174;
t170 = cos(t173);
t185 = g(3) * t170;
t175 = sin(qJ(6));
t177 = sin(qJ(1));
t184 = t177 * t175;
t178 = cos(qJ(6));
t183 = t177 * t178;
t180 = cos(qJ(1));
t182 = t180 * t175;
t181 = t180 * t178;
t167 = g(1) * t177 - g(2) * t180;
t179 = cos(qJ(3));
t176 = sin(qJ(3));
t172 = cos(t174);
t171 = sin(t174);
t169 = sin(t173);
t168 = g(1) * t180 + g(2) * t177;
t1 = [0, -t168, t167, t168, -t167, -g(1) * (t180 * pkin(1) + t177 * qJ(2)) - g(2) * (t177 * pkin(1) - t180 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t179 - t167 * t176, g(3) * t176 - t167 * t179, 0, 0, 0, 0, 0, -g(3) * t172 - t167 * t171, g(3) * t171 - t167 * t172, 0, 0, 0, 0, 0, -t167 * t169 - t185, g(3) * t169 - t167 * t170, 0, 0, 0, 0, 0, -g(1) * (t169 * t183 + t182) - g(2) * (-t169 * t181 + t184) - t178 * t185, -g(1) * (-t169 * t184 + t181) - g(2) * (t169 * t182 + t183) + t175 * t185;];
U_reg  = t1;
