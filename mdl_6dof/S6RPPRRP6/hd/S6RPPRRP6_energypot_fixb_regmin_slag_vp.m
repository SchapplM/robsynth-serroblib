% Calculate minimal parameter regressor of potential energy for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:17
% EndTime: 2019-03-09 02:11:17
% DurationCPUTime: 0.06s
% Computational Cost: add. (65->43), mult. (124->52), div. (0->0), fcn. (127->6), ass. (0->27)
t162 = pkin(2) + pkin(6);
t146 = sin(qJ(4));
t161 = pkin(4) * t146;
t149 = cos(qJ(4));
t160 = g(3) * t149;
t145 = sin(qJ(5));
t147 = sin(qJ(1));
t159 = t147 * t145;
t148 = cos(qJ(5));
t158 = t147 * t148;
t150 = cos(qJ(1));
t157 = t150 * t145;
t156 = t150 * t148;
t155 = t150 * pkin(1) + t147 * qJ(2);
t154 = t150 * qJ(3) + t155;
t153 = t147 * pkin(1) - t150 * qJ(2);
t137 = g(1) * t150 + g(2) * t147;
t152 = t147 * qJ(3) + t153;
t132 = t146 * t159 - t156;
t134 = t146 * t157 + t158;
t151 = g(1) * t134 + g(2) * t132 + t145 * t160;
t136 = g(1) * t147 - g(2) * t150;
t135 = t146 * t156 - t159;
t133 = t146 * t158 + t157;
t131 = -g(3) * t146 + t137 * t149;
t130 = -g(1) * t135 - g(2) * t133 - t148 * t160;
t1 = [0, -t137, t136, t137, -t136, -g(3) * pkin(6) - g(1) * t155 - g(2) * t153, -t136, -t137, -g(1) * t154 - g(2) * t152 - g(3) * t162, 0, 0, 0, 0, 0, -t137 * t146 - t160, -t131, 0, 0, 0, 0, 0, t130, t151, t130, t131, -t151, -g(1) * (t135 * pkin(5) - t147 * pkin(7) + t134 * qJ(6) + t150 * t161 + t154) - g(2) * (t133 * pkin(5) + t150 * pkin(7) + t132 * qJ(6) + t147 * t161 + t152) - g(3) * (t146 * pkin(8) + pkin(3) + t162) + (-g(3) * (pkin(5) * t148 + qJ(6) * t145 + pkin(4)) + t137 * pkin(8)) * t149;];
U_reg  = t1;
