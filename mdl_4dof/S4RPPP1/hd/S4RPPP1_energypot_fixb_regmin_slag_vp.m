% Calculate minimal parameter regressor of potential energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RPPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:35
% EndTime: 2018-11-14 13:45:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (174->44), mult. (197->54), div. (0->0), fcn. (177->10), ass. (0->34)
t154 = pkin(4) + pkin(6);
t160 = sin(t154) / 0.2e1;
t139 = cos(pkin(4));
t159 = t139 * qJ(2) + pkin(5);
t137 = sin(pkin(4));
t140 = sin(qJ(1));
t158 = t137 * t140;
t141 = cos(qJ(1));
t157 = t137 * t141;
t156 = t141 * pkin(1) + qJ(2) * t158;
t155 = pkin(4) - pkin(6);
t150 = sin(t155);
t153 = t160 - t150 / 0.2e1;
t152 = qJ(2) * t157;
t151 = cos(t154);
t149 = g(1) * t140 - g(2) * t141;
t136 = sin(pkin(6));
t147 = cos(t155) / 0.2e1;
t144 = t147 + t151 / 0.2e1;
t118 = t140 * t136 - t141 * t144;
t138 = cos(pkin(6));
t119 = t140 * t138 + t141 * t153;
t134 = t140 * pkin(1);
t148 = t119 * pkin(2) + t118 * qJ(3) + t134;
t126 = t160 + t150 / 0.2e1;
t127 = t147 - t151 / 0.2e1;
t146 = t127 * pkin(2) - t126 * qJ(3) + t159;
t120 = t141 * t136 + t140 * t144;
t121 = t141 * t138 - t140 * t153;
t145 = t121 * pkin(2) + t120 * qJ(3) + t156;
t143 = g(1) * t120 + g(2) * t118 - g(3) * t126;
t142 = g(1) * t121 + g(2) * t119 + g(3) * t127;
t122 = -g(3) * t139 - t149 * t137;
t1 = [0, -g(1) * t141 - g(2) * t140, t149, -t142, t143, t122, -g(1) * t156 - g(2) * (t134 - t152) - g(3) * t159, t122, t142, -t143, -g(1) * t145 - g(2) * (t148 - t152) - g(3) * t146, t122, -t143, -t142, -g(1) * (pkin(3) * t158 + t121 * qJ(4) + t145) - g(2) * (t119 * qJ(4) + (-pkin(3) - qJ(2)) * t157 + t148) - g(3) * (t139 * pkin(3) + t127 * qJ(4) + t146);];
U_reg  = t1;
