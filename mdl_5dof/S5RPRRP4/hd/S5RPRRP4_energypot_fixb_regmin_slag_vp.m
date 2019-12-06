% Calculate minimal parameter regressor of potential energy for
% S5RPRRP4
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:06:51
% EndTime: 2019-12-05 18:06:52
% DurationCPUTime: 0.09s
% Computational Cost: add. (70->41), mult. (97->59), div. (0->0), fcn. (99->8), ass. (0->28)
t160 = qJ(3) + qJ(4);
t155 = cos(t160);
t165 = cos(qJ(3));
t151 = t165 * pkin(3) + pkin(4) * t155 + pkin(2);
t159 = -qJ(5) - pkin(7) - pkin(6);
t161 = sin(pkin(8));
t162 = cos(pkin(8));
t180 = t151 * t162 - t159 * t161;
t179 = g(1) * t161;
t154 = sin(t160);
t164 = sin(qJ(1));
t176 = t164 * t154;
t175 = t164 * t155;
t163 = sin(qJ(3));
t174 = t164 * t163;
t173 = t164 * t165;
t166 = cos(qJ(1));
t172 = t166 * t154;
t171 = t166 * t155;
t170 = t166 * t163;
t169 = t166 * t165;
t168 = t166 * pkin(1) + t164 * qJ(2);
t167 = g(2) * t164 - g(3) * t166;
t157 = t166 * qJ(2);
t153 = g(2) * t166 + g(3) * t164;
t152 = t163 * pkin(3) + pkin(4) * t154;
t150 = g(1) * t162 + t167 * t161;
t1 = [0, t167, t153, t167 * t162 - t179, -t150, -t153, -g(1) * pkin(5) - g(2) * (-t164 * pkin(1) + t157) - g(3) * t168, 0, 0, 0, 0, 0, -t165 * t179 - g(2) * (-t162 * t173 + t170) - g(3) * (t162 * t169 + t174), t163 * t179 - g(2) * (t162 * t174 + t169) - g(3) * (-t162 * t170 + t173), 0, 0, 0, 0, 0, -t155 * t179 - g(2) * (-t162 * t175 + t172) - g(3) * (t162 * t171 + t176), t154 * t179 - g(2) * (t162 * t176 + t171) - g(3) * (-t162 * t172 + t175), t150, -g(1) * (t161 * t151 + t162 * t159 + pkin(5)) - g(2) * (t166 * t152 + t157) - g(3) * (t180 * t166 + t168) + (-g(2) * (-pkin(1) - t180) - g(3) * t152) * t164;];
U_reg = t1;
