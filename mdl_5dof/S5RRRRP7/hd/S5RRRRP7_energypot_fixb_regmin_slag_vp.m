% Calculate minimal parameter regressor of potential energy for
% S5RRRRP7
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
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:43
% EndTime: 2019-12-31 21:57:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (86->35), mult. (109->46), div. (0->0), fcn. (118->8), ass. (0->25)
t164 = qJ(2) + qJ(3);
t163 = cos(t164);
t169 = cos(qJ(2));
t180 = t169 * pkin(2) + pkin(3) * t163 + pkin(1);
t162 = sin(t164);
t178 = g(3) * t162;
t165 = sin(qJ(4));
t167 = sin(qJ(1));
t177 = t167 * t165;
t168 = cos(qJ(4));
t176 = t167 * t168;
t170 = cos(qJ(1));
t175 = t170 * t165;
t174 = t170 * t168;
t173 = g(1) * t170 + g(2) * t167;
t156 = t163 * t177 + t174;
t158 = t163 * t175 - t176;
t172 = g(1) * t158 + g(2) * t156 + t165 * t178;
t171 = -pkin(7) - pkin(6);
t166 = sin(qJ(2));
t159 = t163 * t174 + t177;
t157 = t163 * t176 - t175;
t155 = -g(3) * t163 + t173 * t162;
t154 = -g(1) * t159 - g(2) * t157 - t168 * t178;
t1 = [0, -t173, g(1) * t167 - g(2) * t170, 0, 0, 0, 0, 0, -g(3) * t166 - t173 * t169, -g(3) * t169 + t173 * t166, 0, 0, 0, 0, 0, -t173 * t163 - t178, t155, 0, 0, 0, 0, 0, t154, t172, t154, -t155, -t172, -g(1) * (t159 * pkin(4) + t158 * qJ(5) - t167 * t171 + t180 * t170) - g(2) * (t157 * pkin(4) + t156 * qJ(5) + t180 * t167 + t170 * t171) - g(3) * (t166 * pkin(2) - t163 * pkin(8) + pkin(5)) + (-g(3) * (pkin(4) * t168 + qJ(5) * t165 + pkin(3)) - t173 * pkin(8)) * t162;];
U_reg = t1;
