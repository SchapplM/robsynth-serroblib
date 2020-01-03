% Calculate minimal parameter regressor of potential energy for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:51
% EndTime: 2019-12-31 22:29:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (60->31), mult. (74->51), div. (0->0), fcn. (86->10), ass. (0->22)
t185 = sin(qJ(2));
t198 = g(3) * t185;
t186 = sin(qJ(1));
t188 = cos(qJ(2));
t197 = t186 * t188;
t183 = qJ(3) + qJ(4);
t182 = qJ(5) + t183;
t178 = sin(t182);
t189 = cos(qJ(1));
t196 = t189 * t178;
t179 = cos(t182);
t195 = t189 * t179;
t180 = sin(t183);
t194 = t189 * t180;
t181 = cos(t183);
t193 = t189 * t181;
t184 = sin(qJ(3));
t192 = t189 * t184;
t187 = cos(qJ(3));
t191 = t189 * t187;
t190 = g(1) * t189 + g(2) * t186;
t1 = [0, -t190, g(1) * t186 - g(2) * t189, 0, 0, 0, 0, 0, -t190 * t188 - t198, -g(3) * t188 + t190 * t185, 0, 0, 0, 0, 0, -g(1) * (t186 * t184 + t188 * t191) - g(2) * (t187 * t197 - t192) - t187 * t198, -g(1) * (t186 * t187 - t188 * t192) - g(2) * (-t184 * t197 - t191) + t184 * t198, 0, 0, 0, 0, 0, -g(1) * (t186 * t180 + t188 * t193) - g(2) * (t181 * t197 - t194) - t181 * t198, -g(1) * (t186 * t181 - t188 * t194) - g(2) * (-t180 * t197 - t193) + t180 * t198, 0, 0, 0, 0, 0, -g(1) * (t186 * t178 + t188 * t195) - g(2) * (t179 * t197 - t196) - t179 * t198, -g(1) * (t186 * t179 - t188 * t196) - g(2) * (-t178 * t197 - t195) + t178 * t198;];
U_reg = t1;
