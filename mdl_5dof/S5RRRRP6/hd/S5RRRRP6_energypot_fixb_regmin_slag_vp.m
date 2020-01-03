% Calculate minimal parameter regressor of potential energy for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:33
% EndTime: 2019-12-31 21:54:33
% DurationCPUTime: 0.06s
% Computational Cost: add. (64->31), mult. (76->41), div. (0->0), fcn. (77->8), ass. (0->21)
t151 = qJ(2) + qJ(3);
t149 = sin(t151);
t167 = g(3) * t149;
t153 = sin(qJ(4));
t155 = sin(qJ(1));
t166 = t155 * t153;
t156 = cos(qJ(4));
t165 = t155 * t156;
t158 = cos(qJ(1));
t164 = t158 * t153;
t163 = t158 * t156;
t162 = pkin(4) * t153 + pkin(6) + pkin(7);
t161 = g(1) * t158 + g(2) * t155;
t147 = t156 * pkin(4) + pkin(3);
t150 = cos(t151);
t152 = -qJ(5) - pkin(8);
t157 = cos(qJ(2));
t160 = t157 * pkin(2) + t147 * t150 - t149 * t152 + pkin(1);
t154 = sin(qJ(2));
t146 = -g(3) * t150 + t161 * t149;
t1 = [0, -t161, g(1) * t155 - g(2) * t158, 0, 0, 0, 0, 0, -g(3) * t154 - t161 * t157, -g(3) * t157 + t161 * t154, 0, 0, 0, 0, 0, -t161 * t150 - t167, t146, 0, 0, 0, 0, 0, -g(1) * (t150 * t163 + t166) - g(2) * (t150 * t165 - t164) - t156 * t167, -g(1) * (-t150 * t164 + t165) - g(2) * (-t150 * t166 - t163) + t153 * t167, -t146, -g(3) * (t154 * pkin(2) + t149 * t147 + t150 * t152 + pkin(5)) + (-g(1) * t160 + g(2) * t162) * t158 + (-g(1) * t162 - g(2) * t160) * t155;];
U_reg = t1;
