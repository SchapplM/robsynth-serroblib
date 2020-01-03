% Calculate minimal parameter regressor of potential energy for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:40
% EndTime: 2019-12-31 19:44:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (74->43), mult. (164->63), div. (0->0), fcn. (184->8), ass. (0->25)
t177 = sin(qJ(2));
t192 = qJ(3) * t177 + pkin(1);
t191 = g(3) * t177;
t178 = sin(qJ(1));
t180 = cos(qJ(2));
t189 = t178 * t180;
t174 = sin(pkin(8));
t181 = cos(qJ(1));
t188 = t181 * t174;
t175 = cos(pkin(8));
t187 = t181 * t175;
t186 = t178 * pkin(6) + (pkin(2) * t180 + t192) * t181;
t185 = t177 * pkin(2) - t180 * qJ(3) + pkin(5);
t184 = g(1) * t181 + g(2) * t178;
t183 = pkin(2) * t189 - t181 * pkin(6) + t192 * t178;
t160 = t174 * t189 + t187;
t162 = -t178 * t175 + t180 * t188;
t182 = g(1) * t162 + g(2) * t160 + t174 * t191;
t179 = cos(qJ(5));
t176 = sin(qJ(5));
t163 = t178 * t174 + t180 * t187;
t161 = t175 * t189 - t188;
t159 = -g(3) * t180 + t184 * t177;
t158 = -g(1) * t163 - g(2) * t161 - t175 * t191;
t1 = [0, -t184, g(1) * t178 - g(2) * t181, 0, 0, 0, 0, 0, -t184 * t180 - t191, t159, t158, t182, -t159, -g(1) * t186 - g(2) * t183 - g(3) * t185, t158, -t159, -t182, -g(1) * (t163 * pkin(3) + t162 * qJ(4) + t186) - g(2) * (t161 * pkin(3) + t160 * qJ(4) + t183) - g(3) * ((pkin(3) * t175 + qJ(4) * t174) * t177 + t185), 0, 0, 0, 0, 0, -g(1) * (t162 * t176 + t163 * t179) - g(2) * (t160 * t176 + t161 * t179) - (t174 * t176 + t175 * t179) * t191, -g(1) * (t162 * t179 - t163 * t176) - g(2) * (t160 * t179 - t161 * t176) - (t174 * t179 - t175 * t176) * t191;];
U_reg = t1;
