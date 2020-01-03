% Calculate minimal parameter regressor of potential energy for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:14
% EndTime: 2019-12-31 19:39:14
% DurationCPUTime: 0.09s
% Computational Cost: add. (69->33), mult. (121->46), div. (0->0), fcn. (127->8), ass. (0->24)
t175 = sin(qJ(2));
t192 = qJ(3) * t175 + pkin(1);
t176 = sin(qJ(1));
t178 = cos(qJ(1));
t183 = g(1) * t178 + g(2) * t176;
t177 = cos(qJ(2));
t188 = t176 * t177;
t187 = t177 * t178;
t186 = pkin(2) * t188 + t192 * t176;
t185 = pkin(2) * t187 + t176 * pkin(6) + t192 * t178;
t184 = t175 * pkin(2) - t177 * qJ(3) + pkin(5);
t172 = pkin(8) + qJ(5);
t166 = sin(t172);
t167 = cos(t172);
t182 = t177 * t166 - t175 * t167;
t181 = t175 * t166 + t177 * t167;
t173 = sin(pkin(8));
t174 = cos(pkin(8));
t180 = t177 * t173 - t175 * t174;
t179 = t175 * t173 + t177 * t174;
t161 = g(1) * t176 - g(2) * t178;
t160 = -g(3) * t175 - t183 * t177;
t159 = -g(3) * t177 + t183 * t175;
t1 = [0, -t183, t161, 0, 0, 0, 0, 0, t160, t159, t160, -t161, -t159, -g(1) * t185 - g(2) * (-t178 * pkin(6) + t186) - g(3) * t184, g(3) * t180 - t183 * t179, g(3) * t179 + t183 * t180, t161, -g(1) * (pkin(3) * t187 - t176 * qJ(4) + t185) - g(2) * (pkin(3) * t188 + (-pkin(6) + qJ(4)) * t178 + t186) - g(3) * (t175 * pkin(3) + t184), 0, 0, 0, 0, 0, g(3) * t182 - t183 * t181, g(3) * t181 + t183 * t182;];
U_reg = t1;
