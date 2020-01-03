% Calculate minimal parameter regressor of potential energy for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:01:59
% EndTime: 2019-12-31 21:01:59
% DurationCPUTime: 0.10s
% Computational Cost: add. (92->45), mult. (132->62), div. (0->0), fcn. (138->8), ass. (0->29)
t179 = -qJ(4) - pkin(7);
t181 = sin(qJ(2));
t198 = -t179 * t181 + pkin(1);
t197 = g(3) * t181;
t180 = sin(qJ(3));
t182 = sin(qJ(1));
t195 = t182 * t180;
t184 = cos(qJ(2));
t194 = t182 * t184;
t178 = qJ(3) + pkin(8);
t172 = sin(t178);
t185 = cos(qJ(1));
t193 = t185 * t172;
t173 = cos(t178);
t192 = t185 * t173;
t191 = t185 * t180;
t183 = cos(qJ(3));
t190 = t185 * t183;
t171 = t183 * pkin(3) + pkin(2);
t189 = t181 * t171 + t184 * t179 + pkin(5);
t188 = g(1) * t185 + g(2) * t182;
t187 = pkin(3) * t195 + t182 * pkin(6) + (t171 * t184 + t198) * t185;
t186 = t171 * t194 + (-pkin(3) * t180 - pkin(6)) * t185 + t198 * t182;
t165 = -g(3) * t184 + t188 * t181;
t164 = t182 * t172 + t184 * t192;
t163 = -t182 * t173 + t184 * t193;
t162 = t173 * t194 - t193;
t161 = t172 * t194 + t192;
t1 = [0, -t188, g(1) * t182 - g(2) * t185, 0, 0, 0, 0, 0, -t188 * t184 - t197, t165, 0, 0, 0, 0, 0, -g(1) * (t184 * t190 + t195) - g(2) * (t183 * t194 - t191) - t183 * t197, -g(1) * (t182 * t183 - t184 * t191) - g(2) * (-t180 * t194 - t190) + t180 * t197, -t165, -g(1) * t187 - g(2) * t186 - g(3) * t189, -g(1) * t164 - g(2) * t162 - t173 * t197, -t165, -g(1) * t163 - g(2) * t161 - t172 * t197, -g(1) * (t164 * pkin(4) + t163 * qJ(5) + t187) - g(2) * (t162 * pkin(4) + t161 * qJ(5) + t186) - g(3) * ((pkin(4) * t173 + qJ(5) * t172) * t181 + t189);];
U_reg = t1;
