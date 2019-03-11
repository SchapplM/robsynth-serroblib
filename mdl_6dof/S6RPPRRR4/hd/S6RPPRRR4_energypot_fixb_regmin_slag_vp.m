% Calculate minimal parameter regressor of potential energy for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x30]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:46
% EndTime: 2019-03-09 02:26:46
% DurationCPUTime: 0.07s
% Computational Cost: add. (71->35), mult. (126->60), div. (0->0), fcn. (152->10), ass. (0->24)
t192 = cos(qJ(1));
t191 = sin(qJ(1));
t178 = sin(qJ(4));
t190 = g(3) * t178;
t189 = cos(pkin(10));
t188 = sin(pkin(10));
t176 = qJ(5) + qJ(6);
t171 = sin(t176);
t180 = cos(qJ(4));
t187 = t171 * t180;
t172 = cos(t176);
t186 = t172 * t180;
t177 = sin(qJ(5));
t185 = t177 * t180;
t179 = cos(qJ(5));
t184 = t179 * t180;
t183 = t192 * pkin(1) + t191 * qJ(2);
t164 = -t191 * t188 - t192 * t189;
t165 = t192 * t188 - t191 * t189;
t182 = g(1) * t164 + g(2) * t165;
t181 = t191 * pkin(1) - t192 * qJ(2);
t167 = -g(1) * t192 - g(2) * t191;
t166 = g(1) * t191 - g(2) * t192;
t1 = [0, t167, t166, t167, -t166, -g(3) * pkin(6) - g(1) * t183 - g(2) * t181, t182, g(1) * t165 - g(2) * t164, -g(1) * (t192 * pkin(2) + t183) - g(2) * (t191 * pkin(2) + t181) - g(3) * (-qJ(3) + pkin(6)) 0, 0, 0, 0, 0, t182 * t180 + t190, g(3) * t180 - t182 * t178, 0, 0, 0, 0, 0, -g(1) * (-t164 * t184 + t165 * t177) - g(2) * (-t164 * t177 - t165 * t184) + t179 * t190, -g(1) * (t164 * t185 + t165 * t179) - g(2) * (-t164 * t179 + t165 * t185) - t177 * t190, 0, 0, 0, 0, 0, -g(1) * (-t164 * t186 + t165 * t171) - g(2) * (-t164 * t171 - t165 * t186) + t172 * t190, -g(1) * (t164 * t187 + t165 * t172) - g(2) * (-t164 * t172 + t165 * t187) - t171 * t190;];
U_reg  = t1;
