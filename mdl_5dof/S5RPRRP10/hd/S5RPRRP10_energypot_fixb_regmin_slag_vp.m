% Calculate minimal parameter regressor of potential energy for
% S5RPRRP10
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
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:03
% EndTime: 2019-12-31 18:52:03
% DurationCPUTime: 0.06s
% Computational Cost: add. (69->36), mult. (85->48), div. (0->0), fcn. (83->8), ass. (0->22)
t147 = pkin(8) + qJ(3);
t145 = sin(t147);
t163 = g(3) * t145;
t152 = sin(qJ(4));
t153 = sin(qJ(1));
t162 = t153 * t152;
t154 = cos(qJ(4));
t161 = t153 * t154;
t155 = cos(qJ(1));
t160 = t155 * t152;
t159 = t155 * t154;
t158 = pkin(4) * t152 + pkin(6) + qJ(2);
t157 = g(1) * t155 + g(2) * t153;
t144 = t154 * pkin(4) + pkin(3);
t146 = cos(t147);
t149 = cos(pkin(8));
t150 = -qJ(5) - pkin(7);
t156 = t149 * pkin(2) + t144 * t146 - t145 * t150 + pkin(1);
t148 = sin(pkin(8));
t142 = g(1) * t153 - g(2) * t155;
t141 = -g(3) * t146 + t157 * t145;
t1 = [0, -t157, t142, -g(3) * t148 - t157 * t149, -g(3) * t149 + t157 * t148, -t142, -g(1) * (t155 * pkin(1) + t153 * qJ(2)) - g(2) * (t153 * pkin(1) - t155 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t157 * t146 - t163, t141, 0, 0, 0, 0, 0, -g(1) * (t146 * t159 + t162) - g(2) * (t146 * t161 - t160) - t154 * t163, -g(1) * (-t146 * t160 + t161) - g(2) * (-t146 * t162 - t159) + t152 * t163, -t141, -g(3) * (t148 * pkin(2) + t145 * t144 + t146 * t150 + pkin(5)) + (-g(1) * t156 + g(2) * t158) * t155 + (-g(1) * t158 - g(2) * t156) * t153;];
U_reg = t1;
