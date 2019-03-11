% Calculate minimal parameter regressor of potential energy for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPPRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:41
% EndTime: 2019-03-09 01:37:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (57->37), mult. (105->49), div. (0->0), fcn. (116->8), ass. (0->23)
t146 = pkin(2) + pkin(6);
t133 = sin(qJ(5));
t145 = g(3) * t133;
t144 = sin(pkin(9));
t132 = sin(qJ(6));
t136 = cos(qJ(5));
t143 = t132 * t136;
t135 = cos(qJ(6));
t142 = t135 * t136;
t134 = sin(qJ(1));
t137 = cos(qJ(1));
t141 = t137 * pkin(1) + t134 * qJ(2);
t140 = t137 * qJ(3) + t141;
t128 = t134 * pkin(1);
t139 = -t137 * qJ(2) + t128;
t131 = cos(pkin(9));
t120 = t137 * t131 - t134 * t144;
t121 = t134 * t131 + t137 * t144;
t138 = g(1) * t121 - g(2) * t120;
t125 = t134 * qJ(3);
t123 = g(1) * t137 + g(2) * t134;
t122 = g(1) * t134 - g(2) * t137;
t1 = [0, -t123, t122, t123, -t122, -g(3) * pkin(6) - g(1) * t141 - g(2) * t139, -t122, -t123, -g(1) * t140 - g(2) * (t125 + t139) - g(3) * t146, -t138, -g(1) * t120 - g(2) * t121, -g(1) * (t134 * pkin(3) + t140) - g(2) * (t125 + t128 + (-pkin(3) - qJ(2)) * t137) - g(3) * (qJ(4) + t146) 0, 0, 0, 0, 0, -t138 * t136 - t145, -g(3) * t136 + t138 * t133, 0, 0, 0, 0, 0, -g(1) * (-t120 * t132 + t121 * t142) - g(2) * (-t120 * t142 - t121 * t132) - t135 * t145, -g(1) * (-t120 * t135 - t121 * t143) - g(2) * (t120 * t143 - t121 * t135) + t132 * t145;];
U_reg  = t1;
