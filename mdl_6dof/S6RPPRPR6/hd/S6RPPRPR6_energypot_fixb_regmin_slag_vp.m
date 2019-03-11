% Calculate minimal parameter regressor of potential energy for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:32
% EndTime: 2019-03-09 01:51:32
% DurationCPUTime: 0.06s
% Computational Cost: add. (50->40), mult. (91->48), div. (0->0), fcn. (86->6), ass. (0->22)
t150 = pkin(2) + pkin(6);
t137 = sin(qJ(4));
t149 = pkin(4) * t137;
t148 = g(3) * t137;
t138 = sin(qJ(1));
t140 = cos(qJ(4));
t147 = t138 * t140;
t136 = sin(qJ(6));
t141 = cos(qJ(1));
t146 = t141 * t136;
t139 = cos(qJ(6));
t145 = t141 * t139;
t144 = t141 * pkin(1) + t138 * qJ(2);
t143 = t141 * qJ(3) + t144;
t133 = t138 * pkin(1);
t142 = -t141 * qJ(2) + t133;
t129 = g(1) * t141 + g(2) * t138;
t130 = t138 * qJ(3);
t128 = g(1) * t138 - g(2) * t141;
t127 = t129 * t140 - t148;
t126 = g(3) * t140 + t129 * t137;
t1 = [0, -t129, t128, t129, -t128, -g(3) * pkin(6) - g(1) * t144 - g(2) * t142, -t128, -t129, -g(1) * t143 - g(2) * (t130 + t142) - g(3) * t150, 0, 0, 0, 0, 0, -t126, -t127, t128, t126, t127, -g(1) * (-t138 * pkin(7) + t143) - g(2) * (-qJ(5) * t147 + t138 * t149 + t130 + t133) - g(3) * (t140 * pkin(4) + t137 * qJ(5) + pkin(3) + t150) + (-g(1) * (-qJ(5) * t140 + t149) - g(2) * (pkin(7) - qJ(2))) * t141, 0, 0, 0, 0, 0, -g(1) * (-t138 * t139 - t140 * t146) - g(2) * (-t136 * t147 + t145) - t136 * t148, -g(1) * (t138 * t136 - t140 * t145) - g(2) * (-t139 * t147 - t146) - t139 * t148;];
U_reg  = t1;
