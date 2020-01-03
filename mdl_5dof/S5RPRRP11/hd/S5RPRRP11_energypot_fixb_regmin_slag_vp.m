% Calculate minimal parameter regressor of potential energy for
% S5RPRRP11
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
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP11_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP11_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:54:31
% EndTime: 2019-12-31 18:54:31
% DurationCPUTime: 0.08s
% Computational Cost: add. (91->40), mult. (118->53), div. (0->0), fcn. (124->8), ass. (0->26)
t160 = pkin(8) + qJ(3);
t159 = cos(t160);
t162 = cos(pkin(8));
t176 = t162 * pkin(2) + pkin(3) * t159 + pkin(1);
t158 = sin(t160);
t174 = g(3) * t158;
t164 = sin(qJ(4));
t165 = sin(qJ(1));
t173 = t165 * t164;
t166 = cos(qJ(4));
t172 = t165 * t166;
t167 = cos(qJ(1));
t171 = t167 * t164;
t170 = t167 * t166;
t169 = g(1) * t167 + g(2) * t165;
t151 = t159 * t173 + t170;
t153 = t159 * t171 - t172;
t168 = g(1) * t153 + g(2) * t151 + t164 * t174;
t163 = -pkin(6) - qJ(2);
t161 = sin(pkin(8));
t155 = g(1) * t165 - g(2) * t167;
t154 = t159 * t170 + t173;
t152 = t159 * t172 - t171;
t150 = -g(3) * t159 + t169 * t158;
t149 = -g(1) * t154 - g(2) * t152 - t166 * t174;
t1 = [0, -t169, t155, -g(3) * t161 - t169 * t162, -g(3) * t162 + t169 * t161, -t155, -g(1) * (t167 * pkin(1) + t165 * qJ(2)) - g(2) * (t165 * pkin(1) - t167 * qJ(2)) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t169 * t159 - t174, t150, 0, 0, 0, 0, 0, t149, t168, t149, -t150, -t168, -g(1) * (t154 * pkin(4) + t153 * qJ(5) - t165 * t163 + t176 * t167) - g(2) * (t152 * pkin(4) + t151 * qJ(5) + t167 * t163 + t176 * t165) - g(3) * (t161 * pkin(2) - t159 * pkin(7) + pkin(5)) + (-g(3) * (pkin(4) * t166 + qJ(5) * t164 + pkin(3)) - t169 * pkin(7)) * t158;];
U_reg = t1;
