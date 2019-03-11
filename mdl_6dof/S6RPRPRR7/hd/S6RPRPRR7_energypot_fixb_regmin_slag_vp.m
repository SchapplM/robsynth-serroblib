% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:56:03
% EndTime: 2019-03-09 03:56:03
% DurationCPUTime: 0.06s
% Computational Cost: add. (61->32), mult. (71->42), div. (0->0), fcn. (69->8), ass. (0->21)
t182 = sin(qJ(3));
t193 = pkin(3) * t182;
t176 = qJ(3) + pkin(10) + qJ(5);
t175 = cos(t176);
t192 = g(3) * t175;
t181 = sin(qJ(6));
t183 = sin(qJ(1));
t191 = t183 * t181;
t184 = cos(qJ(6));
t190 = t183 * t184;
t186 = cos(qJ(1));
t189 = t186 * t181;
t188 = t186 * t184;
t187 = t186 * pkin(1) + t183 * qJ(2);
t172 = g(1) * t183 - g(2) * t186;
t185 = cos(qJ(3));
t180 = -qJ(4) - pkin(7);
t178 = t183 * pkin(1);
t174 = sin(t176);
t173 = g(1) * t186 + g(2) * t183;
t1 = [0, -t173, t172, t173, -t172, -g(1) * t187 - g(2) * (-t186 * qJ(2) + t178) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t185 - t172 * t182, g(3) * t182 - t172 * t185, -t173, -g(1) * (-t186 * t180 + t183 * t193 + t187) - g(2) * (-t183 * t180 + t178 + (-qJ(2) - t193) * t186) - g(3) * (t185 * pkin(3) + pkin(2) + pkin(6)) 0, 0, 0, 0, 0, -t172 * t174 - t192, g(3) * t174 - t172 * t175, 0, 0, 0, 0, 0, -g(1) * (t174 * t190 + t189) - g(2) * (-t174 * t188 + t191) - t184 * t192, -g(1) * (-t174 * t191 + t188) - g(2) * (t174 * t189 + t190) + t181 * t192;];
U_reg  = t1;
