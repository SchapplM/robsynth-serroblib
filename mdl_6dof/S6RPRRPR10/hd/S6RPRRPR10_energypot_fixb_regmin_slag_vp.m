% Calculate minimal parameter regressor of potential energy for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% U_reg [1x31]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:37:09
% EndTime: 2019-03-09 05:37:09
% DurationCPUTime: 0.08s
% Computational Cost: add. (71->46), mult. (155->65), div. (0->0), fcn. (175->8), ass. (0->26)
t184 = sin(qJ(3));
t198 = pkin(3) * t184;
t188 = cos(qJ(3));
t197 = g(3) * t188;
t183 = sin(qJ(4));
t185 = sin(qJ(1));
t196 = t185 * t183;
t187 = cos(qJ(4));
t195 = t185 * t187;
t189 = cos(qJ(1));
t194 = t189 * t183;
t193 = t189 * t187;
t192 = t189 * pkin(1) + t185 * qJ(2);
t191 = t185 * pkin(1) - t189 * qJ(2);
t176 = g(1) * t185 - g(2) * t189;
t172 = t184 * t196 - t193;
t174 = t184 * t194 + t195;
t190 = g(1) * t172 - g(2) * t174 + t183 * t197;
t186 = cos(qJ(6));
t182 = sin(qJ(6));
t177 = g(1) * t189 + g(2) * t185;
t175 = -t184 * t193 + t196;
t173 = t184 * t195 + t194;
t171 = -g(3) * t184 + t176 * t188;
t170 = -g(1) * t173 - g(2) * t175 - t187 * t197;
t1 = [0, -t177, t176, t177, -t176, -g(3) * pkin(6) - g(1) * t192 - g(2) * t191, 0, 0, 0, 0, 0, -t176 * t184 - t197, -t171, 0, 0, 0, 0, 0, t170, t190, t170, t171, -t190, -g(1) * (t173 * pkin(4) + t189 * pkin(7) + t172 * qJ(5) + t185 * t198 + t192) - g(2) * (t175 * pkin(4) + t185 * pkin(7) - t174 * qJ(5) - t189 * t198 + t191) - g(3) * (t184 * pkin(8) + pkin(2) + pkin(6)) + (-g(3) * (pkin(4) * t187 + qJ(5) * t183 + pkin(3)) + t176 * pkin(8)) * t188, 0, 0, 0, 0, 0, -g(1) * (t172 * t182 + t173 * t186) - g(2) * (-t174 * t182 + t175 * t186) - (t182 * t183 + t186 * t187) * t197, -g(1) * (t172 * t186 - t173 * t182) - g(2) * (-t174 * t186 - t175 * t182) - (-t182 * t187 + t183 * t186) * t197;];
U_reg  = t1;
