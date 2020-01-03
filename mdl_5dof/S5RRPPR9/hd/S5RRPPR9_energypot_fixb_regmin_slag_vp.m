% Calculate minimal parameter regressor of potential energy for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:41:47
% EndTime: 2019-12-31 19:41:47
% DurationCPUTime: 0.07s
% Computational Cost: add. (49->35), mult. (99->44), div. (0->0), fcn. (97->6), ass. (0->22)
t141 = sin(qJ(2));
t158 = qJ(3) * t141 + pkin(1);
t144 = cos(qJ(2));
t157 = g(3) * t144;
t140 = sin(qJ(5));
t142 = sin(qJ(1));
t155 = t142 * t140;
t143 = cos(qJ(5));
t154 = t142 * t143;
t153 = t142 * t144;
t145 = cos(qJ(1));
t152 = t144 * t145;
t151 = t145 * t140;
t150 = t145 * t143;
t149 = pkin(2) * t153 + t158 * t142;
t148 = pkin(2) * t152 + t142 * pkin(6) + t158 * t145;
t147 = t141 * pkin(2) - t144 * qJ(3) + pkin(5);
t146 = g(1) * t145 + g(2) * t142;
t130 = g(1) * t142 - g(2) * t145;
t129 = g(3) * t141 + t146 * t144;
t128 = t146 * t141 - t157;
t1 = [0, -t146, t130, 0, 0, 0, 0, 0, -t129, t128, -t129, -t130, -t128, -g(1) * t148 - g(2) * (-t145 * pkin(6) + t149) - g(3) * t147, -t128, t129, t130, -g(1) * (pkin(3) * t152 - t142 * qJ(4) + t148) - g(2) * (pkin(3) * t153 + (-pkin(6) + qJ(4)) * t145 + t149) - g(3) * (t141 * pkin(3) + t147), 0, 0, 0, 0, 0, -g(1) * (t141 * t150 - t155) - g(2) * (t141 * t154 + t151) + t143 * t157, -g(1) * (-t141 * t151 - t154) - g(2) * (-t141 * t155 + t150) - t140 * t157;];
U_reg = t1;
