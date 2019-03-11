% Calculate minimal parameter regressor of potential energy for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:59:59
% EndTime: 2019-03-09 02:59:59
% DurationCPUTime: 0.06s
% Computational Cost: add. (61->43), mult. (114->55), div. (0->0), fcn. (109->6), ass. (0->23)
t147 = cos(qJ(1));
t158 = g(2) * t147;
t143 = sin(qJ(3));
t157 = g(3) * t143;
t144 = sin(qJ(1));
t156 = t143 * t144;
t146 = cos(qJ(3));
t155 = t144 * t146;
t142 = sin(qJ(6));
t154 = t147 * t142;
t145 = cos(qJ(6));
t153 = t147 * t145;
t152 = t147 * pkin(1) + t144 * qJ(2);
t138 = t144 * pkin(1);
t151 = t147 * t146 * qJ(4) + t144 * pkin(7) + t138;
t150 = -pkin(3) * t143 - qJ(2);
t149 = t146 * pkin(3) + t143 * qJ(4) + pkin(2) + pkin(6);
t130 = g(1) * t144 - t158;
t148 = pkin(3) * t156 + t147 * pkin(7) - qJ(4) * t155 + t152;
t131 = g(1) * t147 + g(2) * t144;
t129 = t130 * t146 - t157;
t128 = g(1) * t156 + g(3) * t146 - t143 * t158;
t1 = [0, -t131, t130, t131, -t130, -g(1) * t152 - g(2) * (-t147 * qJ(2) + t138) - g(3) * pkin(6), 0, 0, 0, 0, 0, -t128, -t129, -t128, -t131, t129, -g(1) * t148 - g(2) * (t150 * t147 + t151) - g(3) * t149, t129, t128, t131, -g(1) * (pkin(4) * t156 + t148) - g(2) * (-t144 * qJ(5) + t151) - g(3) * (t146 * pkin(4) + t149) + (g(1) * qJ(5) - g(2) * (-pkin(4) * t143 + t150)) * t147, 0, 0, 0, 0, 0, -g(1) * (-t145 * t155 - t154) - g(2) * (-t144 * t142 + t146 * t153) - t145 * t157, -g(1) * (t142 * t155 - t153) - g(2) * (-t144 * t145 - t146 * t154) + t142 * t157;];
U_reg  = t1;
