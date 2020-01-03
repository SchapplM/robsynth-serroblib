% Calculate minimal parameter regressor of potential energy for
% S5RRRRP8
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
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:08
% EndTime: 2019-12-31 22:02:08
% DurationCPUTime: 0.07s
% Computational Cost: add. (65->36), mult. (88->52), div. (0->0), fcn. (93->8), ass. (0->22)
t162 = sin(qJ(2));
t175 = g(3) * t162;
t160 = qJ(3) + qJ(4);
t157 = sin(t160);
t161 = sin(qJ(3));
t174 = t161 * pkin(3) + pkin(4) * t157 + pkin(6);
t163 = sin(qJ(1));
t165 = cos(qJ(2));
t173 = t163 * t165;
t166 = cos(qJ(1));
t172 = t166 * t157;
t158 = cos(t160);
t171 = t166 * t158;
t170 = t166 * t161;
t164 = cos(qJ(3));
t169 = t166 * t164;
t168 = g(1) * t166 + g(2) * t163;
t155 = t164 * pkin(3) + pkin(4) * t158 + pkin(2);
t159 = -qJ(5) - pkin(8) - pkin(7);
t167 = t155 * t165 - t159 * t162 + pkin(1);
t154 = -g(3) * t165 + t168 * t162;
t1 = [0, -t168, g(1) * t163 - g(2) * t166, 0, 0, 0, 0, 0, -t168 * t165 - t175, t154, 0, 0, 0, 0, 0, -g(1) * (t163 * t161 + t165 * t169) - g(2) * (t164 * t173 - t170) - t164 * t175, -g(1) * (t163 * t164 - t165 * t170) - g(2) * (-t161 * t173 - t169) + t161 * t175, 0, 0, 0, 0, 0, -g(1) * (t163 * t157 + t165 * t171) - g(2) * (t158 * t173 - t172) - t158 * t175, -g(1) * (t163 * t158 - t165 * t172) - g(2) * (-t157 * t173 - t171) + t157 * t175, -t154, -g(3) * (t162 * t155 + t165 * t159 + pkin(5)) + (-g(1) * t167 + g(2) * t174) * t166 + (-g(1) * t174 - g(2) * t167) * t163;];
U_reg = t1;
