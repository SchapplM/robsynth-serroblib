% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x34]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:21:12
% EndTime: 2019-03-09 07:21:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (60->30), mult. (75->47), div. (0->0), fcn. (80->10), ass. (0->24)
t180 = qJ(3) + qJ(4);
t178 = cos(t180);
t195 = g(3) * t178;
t179 = qJ(5) + qJ(6);
t175 = sin(t179);
t183 = sin(qJ(1));
t194 = t183 * t175;
t177 = cos(t179);
t193 = t183 * t177;
t181 = sin(qJ(5));
t192 = t183 * t181;
t184 = cos(qJ(5));
t191 = t183 * t184;
t186 = cos(qJ(1));
t190 = t186 * t175;
t189 = t186 * t177;
t188 = t186 * t181;
t187 = t186 * t184;
t173 = g(1) * t183 - g(2) * t186;
t185 = cos(qJ(3));
t182 = sin(qJ(3));
t176 = sin(t180);
t174 = g(1) * t186 + g(2) * t183;
t1 = [0, -t174, t173, t174, -t173, -g(1) * (t186 * pkin(1) + t183 * qJ(2)) - g(2) * (t183 * pkin(1) - t186 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t185 - t173 * t182, g(3) * t182 - t173 * t185, 0, 0, 0, 0, 0, -t173 * t176 - t195, g(3) * t176 - t173 * t178, 0, 0, 0, 0, 0, -g(1) * (t176 * t191 + t188) - g(2) * (-t176 * t187 + t192) - t184 * t195, -g(1) * (-t176 * t192 + t187) - g(2) * (t176 * t188 + t191) + t181 * t195, 0, 0, 0, 0, 0, -g(1) * (t176 * t193 + t190) - g(2) * (-t176 * t189 + t194) - t177 * t195, -g(1) * (-t176 * t194 + t189) - g(2) * (t176 * t190 + t193) + t175 * t195;];
U_reg  = t1;
