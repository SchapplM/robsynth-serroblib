% Calculate minimal parameter regressor of potential energy for
% S6RPPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:18
% EndTime: 2019-03-09 02:06:18
% DurationCPUTime: 0.08s
% Computational Cost: add. (108->44), mult. (210->63), div. (0->0), fcn. (255->8), ass. (0->30)
t193 = cos(qJ(1));
t192 = sin(qJ(1));
t176 = sin(qJ(4));
t191 = g(3) * t176;
t190 = -qJ(3) + pkin(6);
t189 = cos(pkin(9));
t188 = sin(pkin(9));
t175 = sin(qJ(5));
t178 = cos(qJ(4));
t187 = t175 * t178;
t177 = cos(qJ(5));
t186 = t177 * t178;
t185 = t193 * pkin(1) + t192 * qJ(2);
t184 = t193 * pkin(2) + t185;
t163 = -t192 * t188 - t193 * t189;
t164 = t193 * t188 - t192 * t189;
t183 = g(1) * t163 + g(2) * t164;
t182 = t192 * pkin(1) - t193 * qJ(2);
t181 = pkin(4) * t178 + pkin(8) * t176 + pkin(3);
t180 = t192 * pkin(2) + t182;
t157 = t163 * t177 - t164 * t187;
t159 = -t163 * t187 - t164 * t177;
t179 = -g(1) * t159 - g(2) * t157 + t175 * t191;
t166 = -g(1) * t193 - g(2) * t192;
t165 = g(1) * t192 - g(2) * t193;
t160 = -t163 * t186 + t164 * t175;
t158 = -t163 * t175 - t164 * t186;
t156 = g(3) * t178 - t183 * t176;
t155 = -g(1) * t160 - g(2) * t158 + t177 * t191;
t1 = [0, t166, t165, t166, -t165, -g(3) * pkin(6) - g(1) * t185 - g(2) * t182, t183, g(1) * t164 - g(2) * t163, -g(1) * t184 - g(2) * t180 - g(3) * t190, 0, 0, 0, 0, 0, t183 * t178 + t191, t156, 0, 0, 0, 0, 0, t155, -t179, t155, -t156, t179, -g(1) * (t160 * pkin(5) + t159 * qJ(6) + t184) - g(2) * (t158 * pkin(5) + t157 * qJ(6) + t180) - g(3) * (t178 * pkin(8) + t190) - (-pkin(5) * t177 - qJ(6) * t175 - pkin(4)) * t191 + (-g(1) * pkin(7) + g(2) * t181) * t164 + (g(2) * pkin(7) + g(1) * t181) * t163;];
U_reg  = t1;
