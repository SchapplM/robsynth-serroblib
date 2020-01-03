% Calculate minimal parameter regressor of potential energy for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:52
% EndTime: 2019-12-31 18:32:52
% DurationCPUTime: 0.06s
% Computational Cost: add. (67->35), mult. (85->46), div. (0->0), fcn. (83->8), ass. (0->21)
t150 = pkin(8) + qJ(3);
t148 = cos(t150);
t164 = g(3) * t148;
t154 = sin(qJ(5));
t155 = sin(qJ(1));
t163 = t155 * t154;
t156 = cos(qJ(5));
t162 = t155 * t156;
t157 = cos(qJ(1));
t161 = t157 * t154;
t160 = t157 * t156;
t159 = g(1) * t157 + g(2) * t155;
t147 = sin(t150);
t152 = cos(pkin(8));
t158 = t152 * pkin(2) + pkin(3) * t148 + qJ(4) * t147 + pkin(1);
t153 = -pkin(6) - qJ(2);
t151 = sin(pkin(8));
t145 = g(1) * t155 - g(2) * t157;
t144 = g(3) * t147 + t159 * t148;
t143 = t159 * t147 - t164;
t1 = [0, -t159, t145, -g(3) * t151 - t159 * t152, -g(3) * t152 + t159 * t151, -t145, -g(1) * (t157 * pkin(1) + t155 * qJ(2)) - g(2) * (t155 * pkin(1) - t157 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t144, t143, -t145, t144, -t143, -g(3) * (t151 * pkin(2) + t147 * pkin(3) - t148 * qJ(4) + pkin(5)) + (-g(1) * t158 - g(2) * t153) * t157 + (g(1) * t153 - g(2) * t158) * t155, 0, 0, 0, 0, 0, -g(1) * (t147 * t161 + t162) - g(2) * (t147 * t163 - t160) + t154 * t164, -g(1) * (t147 * t160 - t163) - g(2) * (t147 * t162 + t161) + t156 * t164;];
U_reg = t1;
