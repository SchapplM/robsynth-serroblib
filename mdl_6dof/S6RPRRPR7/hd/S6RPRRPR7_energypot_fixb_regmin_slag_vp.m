% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:21:00
% EndTime: 2019-03-09 05:21:00
% DurationCPUTime: 0.06s
% Computational Cost: add. (63->35), mult. (74->45), div. (0->0), fcn. (72->10), ass. (0->23)
t183 = qJ(3) + qJ(4);
t176 = pkin(10) + t183;
t195 = g(3) * cos(t176);
t184 = sin(qJ(6));
t186 = sin(qJ(1));
t194 = t186 * t184;
t187 = cos(qJ(6));
t193 = t186 * t187;
t189 = cos(qJ(1));
t192 = t189 * t184;
t191 = t189 * t187;
t190 = t189 * pkin(1) + t186 * qJ(2);
t172 = g(1) * t186 - g(2) * t189;
t188 = cos(qJ(3));
t185 = sin(qJ(3));
t182 = -qJ(5) - pkin(8) - pkin(7);
t180 = t186 * pkin(1);
t178 = cos(t183);
t177 = sin(t183);
t174 = sin(t176);
t173 = g(1) * t189 + g(2) * t186;
t171 = t185 * pkin(3) + pkin(4) * t177;
t1 = [0, -t173, t172, t173, -t172, -g(1) * t190 - g(2) * (-t189 * qJ(2) + t180) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t188 - t172 * t185, g(3) * t185 - t172 * t188, 0, 0, 0, 0, 0, -g(3) * t178 - t172 * t177, g(3) * t177 - t172 * t178, -t173, -g(1) * (t186 * t171 - t189 * t182 + t190) - g(2) * (-t186 * t182 + t180 + (-qJ(2) - t171) * t189) - g(3) * (t188 * pkin(3) + pkin(4) * t178 + pkin(2) + pkin(6)) 0, 0, 0, 0, 0, -g(1) * (t174 * t193 + t192) - g(2) * (-t174 * t191 + t194) - t187 * t195, -g(1) * (-t174 * t194 + t191) - g(2) * (t174 * t192 + t193) + t184 * t195;];
U_reg  = t1;
