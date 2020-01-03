% Calculate minimal parameter regressor of potential energy for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:36
% EndTime: 2019-12-31 18:45:36
% DurationCPUTime: 0.06s
% Computational Cost: add. (88->34), mult. (103->49), div. (0->0), fcn. (109->8), ass. (0->24)
t146 = sin(qJ(3));
t158 = g(3) * t146;
t157 = qJ(2) + pkin(5);
t145 = sin(qJ(4));
t149 = cos(qJ(3));
t156 = t145 * t149;
t148 = cos(qJ(4));
t155 = t148 * t149;
t144 = qJ(1) + pkin(8);
t142 = sin(t144);
t143 = cos(t144);
t154 = g(1) * t143 + g(2) * t142;
t147 = sin(qJ(1));
t150 = cos(qJ(1));
t153 = -g(1) * t150 - g(2) * t147;
t152 = pkin(3) * t149 + pkin(7) * t146 + pkin(2);
t137 = t142 * t156 + t143 * t148;
t139 = -t142 * t148 + t143 * t156;
t151 = g(1) * t139 + g(2) * t137 + t145 * t158;
t140 = t142 * t145 + t143 * t155;
t138 = t142 * t155 - t143 * t145;
t136 = -g(3) * t149 + t154 * t146;
t135 = -g(1) * t140 - g(2) * t138 - t148 * t158;
t1 = [0, t153, g(1) * t147 - g(2) * t150, t153 * pkin(1) - g(3) * t157, 0, 0, 0, 0, 0, -t154 * t149 - t158, t136, 0, 0, 0, 0, 0, t135, t151, t135, -t136, -t151, -g(1) * (t150 * pkin(1) + t140 * pkin(4) + t139 * qJ(5)) - g(2) * (t147 * pkin(1) + t138 * pkin(4) + t137 * qJ(5)) - g(3) * (-t149 * pkin(7) + t157) - (pkin(4) * t148 + qJ(5) * t145 + pkin(3)) * t158 + (g(2) * pkin(6) - g(1) * t152) * t143 + (-g(1) * pkin(6) - g(2) * t152) * t142;];
U_reg = t1;
