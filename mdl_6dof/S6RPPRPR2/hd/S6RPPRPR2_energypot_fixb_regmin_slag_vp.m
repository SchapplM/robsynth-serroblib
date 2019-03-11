% Calculate minimal parameter regressor of potential energy for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:36
% EndTime: 2019-03-09 01:42:36
% DurationCPUTime: 0.07s
% Computational Cost: add. (112->41), mult. (94->55), div. (0->0), fcn. (89->10), ass. (0->29)
t177 = pkin(10) + qJ(4);
t173 = cos(t177);
t195 = g(3) * t173;
t181 = qJ(2) + pkin(6);
t194 = g(3) * t181;
t178 = qJ(1) + pkin(9);
t172 = sin(t178);
t183 = sin(qJ(6));
t193 = t172 * t183;
t185 = cos(qJ(6));
t192 = t172 * t185;
t174 = cos(t178);
t191 = t174 * t183;
t190 = t174 * t185;
t189 = g(1) * t174 + g(2) * t172;
t184 = sin(qJ(1));
t186 = cos(qJ(1));
t188 = -g(1) * t186 - g(2) * t184;
t171 = sin(t177);
t180 = cos(pkin(10));
t187 = t180 * pkin(3) + pkin(4) * t173 + qJ(5) * t171 + pkin(2);
t182 = -pkin(7) - qJ(3);
t179 = sin(pkin(10));
t176 = t186 * pkin(1);
t175 = t184 * pkin(1);
t169 = -g(1) * t172 + g(2) * t174;
t168 = g(3) * t171 + t189 * t173;
t167 = t189 * t171 - t195;
t1 = [0, t188, g(1) * t184 - g(2) * t186, t188 * pkin(1) - t194, -g(3) * t179 - t189 * t180, -g(3) * t180 + t189 * t179, t169, -g(1) * (t174 * pkin(2) + t172 * qJ(3) + t176) - g(2) * (t172 * pkin(2) - t174 * qJ(3) + t175) - t194, 0, 0, 0, 0, 0, -t168, t167, t169, t168, -t167, -g(1) * t176 - g(2) * t175 - g(3) * (t179 * pkin(3) + t171 * pkin(4) - t173 * qJ(5) + t181) + (-g(1) * t187 - g(2) * t182) * t174 + (g(1) * t182 - g(2) * t187) * t172, 0, 0, 0, 0, 0, -g(1) * (t171 * t191 + t192) - g(2) * (t171 * t193 - t190) + t183 * t195, -g(1) * (t171 * t190 - t193) - g(2) * (t171 * t192 + t191) + t185 * t195;];
U_reg  = t1;
