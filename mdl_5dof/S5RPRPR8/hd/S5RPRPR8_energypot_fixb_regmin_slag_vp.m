% Calculate minimal parameter regressor of potential energy for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:10
% EndTime: 2019-12-31 18:22:11
% DurationCPUTime: 0.07s
% Computational Cost: add. (80->35), mult. (83->55), div. (0->0), fcn. (85->10), ass. (0->24)
t155 = sin(qJ(3));
t168 = g(3) * t155;
t167 = qJ(2) + pkin(5);
t152 = qJ(1) + pkin(8);
t148 = sin(t152);
t157 = cos(qJ(3));
t166 = t148 * t157;
t150 = cos(t152);
t165 = t150 * t157;
t153 = sin(pkin(9));
t164 = t153 * t157;
t154 = cos(pkin(9));
t163 = t154 * t157;
t162 = g(1) * t150 + g(2) * t148;
t156 = sin(qJ(1));
t158 = cos(qJ(1));
t161 = -g(1) * t158 - g(2) * t156;
t160 = pkin(3) * t157 + qJ(4) * t155 + pkin(2);
t159 = t161 * pkin(1);
t151 = pkin(9) + qJ(5);
t149 = cos(t151);
t147 = sin(t151);
t146 = -g(3) * t157 + t162 * t155;
t1 = [0, t161, g(1) * t156 - g(2) * t158, -g(3) * t167 + t159, 0, 0, 0, 0, 0, -t162 * t157 - t168, t146, -g(1) * (t148 * t153 + t150 * t163) - g(2) * (t148 * t163 - t150 * t153) - t154 * t168, -g(1) * (t148 * t154 - t150 * t164) - g(2) * (-t148 * t164 - t150 * t154) + t153 * t168, -t146, -g(3) * (t155 * pkin(3) - t157 * qJ(4) + t167) + t159 + (g(2) * pkin(6) - g(1) * t160) * t150 + (-g(1) * pkin(6) - g(2) * t160) * t148, 0, 0, 0, 0, 0, -g(1) * (t148 * t147 + t149 * t165) - g(2) * (-t150 * t147 + t149 * t166) - t149 * t168, -g(1) * (-t147 * t165 + t148 * t149) - g(2) * (-t147 * t166 - t150 * t149) + t147 * t168;];
U_reg = t1;
