% Calculate minimal parameter regressor of potential energy for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:51
% EndTime: 2019-03-09 01:35:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (66->34), mult. (123->50), div. (0->0), fcn. (142->8), ass. (0->23)
t157 = cos(qJ(1));
t156 = sin(qJ(1));
t155 = g(3) * (-qJ(3) + pkin(6));
t143 = cos(qJ(5));
t154 = g(3) * t143;
t153 = cos(pkin(9));
t152 = sin(pkin(9));
t140 = sin(qJ(6));
t141 = sin(qJ(5));
t151 = t140 * t141;
t142 = cos(qJ(6));
t150 = t141 * t142;
t149 = t157 * pkin(1) + t156 * qJ(2);
t148 = t157 * pkin(2) + t149;
t127 = -t156 * t152 - t157 * t153;
t128 = t157 * t152 - t156 * t153;
t147 = g(1) * t128 - g(2) * t127;
t146 = g(1) * t127 + g(2) * t128;
t145 = t156 * pkin(1) - t157 * qJ(2);
t144 = t156 * pkin(2) + t145;
t130 = -g(1) * t157 - g(2) * t156;
t129 = g(1) * t156 - g(2) * t157;
t1 = [0, t130, t129, t130, -t129, -g(3) * pkin(6) - g(1) * t149 - g(2) * t145, t146, t147, -g(1) * t148 - g(2) * t144 - t155, -t146, -t147, -g(1) * (-t127 * pkin(3) + t128 * qJ(4) + t148) - g(2) * (-t128 * pkin(3) - t127 * qJ(4) + t144) - t155, 0, 0, 0, 0, 0, -t147 * t141 + t154, -g(3) * t141 - t147 * t143, 0, 0, 0, 0, 0, -g(1) * (-t127 * t140 + t128 * t150) - g(2) * (-t127 * t150 - t128 * t140) + t142 * t154, -g(1) * (-t127 * t142 - t128 * t151) - g(2) * (t127 * t151 - t128 * t142) - t140 * t154;];
U_reg  = t1;
