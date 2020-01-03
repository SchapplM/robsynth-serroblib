% Calculate minimal parameter regressor of potential energy for
% S5RRRRP9
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:03
% EndTime: 2019-12-31 22:06:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (90->41), mult. (123->58), div. (0->0), fcn. (136->8), ass. (0->28)
t171 = sin(qJ(2));
t186 = g(3) * t171;
t172 = sin(qJ(1));
t174 = cos(qJ(2));
t185 = t172 * t174;
t169 = qJ(3) + qJ(4);
t167 = sin(t169);
t175 = cos(qJ(1));
t184 = t175 * t167;
t168 = cos(t169);
t183 = t175 * t168;
t170 = sin(qJ(3));
t182 = t175 * t170;
t173 = cos(qJ(3));
t181 = t175 * t173;
t180 = pkin(3) * t170 + pkin(6);
t179 = g(1) * t175 + g(2) * t172;
t166 = t173 * pkin(3) + pkin(2);
t176 = -pkin(8) - pkin(7);
t178 = t166 * t174 - t171 * t176 + pkin(1);
t160 = t167 * t185 + t183;
t162 = -t172 * t168 + t174 * t184;
t177 = g(1) * t162 + g(2) * t160 + t167 * t186;
t164 = -g(3) * t174 + t179 * t171;
t163 = t172 * t167 + t174 * t183;
t161 = t168 * t185 - t184;
t159 = -g(1) * t163 - g(2) * t161 - t168 * t186;
t1 = [0, -t179, g(1) * t172 - g(2) * t175, 0, 0, 0, 0, 0, -t179 * t174 - t186, t164, 0, 0, 0, 0, 0, -g(1) * (t172 * t170 + t174 * t181) - g(2) * (t173 * t185 - t182) - t173 * t186, -g(1) * (t172 * t173 - t174 * t182) - g(2) * (-t170 * t185 - t181) + t170 * t186, 0, 0, 0, 0, 0, t159, t177, t159, -t164, -t177, -g(1) * (t163 * pkin(4) + t162 * qJ(5)) - g(2) * (t161 * pkin(4) + t160 * qJ(5)) - g(3) * (t174 * t176 + pkin(5)) - (pkin(4) * t168 + qJ(5) * t167 + t166) * t186 + (-g(1) * t178 + g(2) * t180) * t175 + (-g(1) * t180 - g(2) * t178) * t172;];
U_reg = t1;
