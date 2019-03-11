% Calculate minimal parameter regressor of potential energy for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:13:00
% EndTime: 2019-03-09 03:13:00
% DurationCPUTime: 0.09s
% Computational Cost: add. (129->44), mult. (141->58), div. (0->0), fcn. (144->8), ass. (0->28)
t186 = sin(qJ(3));
t189 = cos(qJ(3));
t203 = pkin(3) * t189 + qJ(4) * t186 + pkin(2);
t201 = g(3) * t189;
t200 = qJ(2) + pkin(6);
t185 = sin(qJ(5));
t198 = t185 * t186;
t188 = cos(qJ(5));
t197 = t186 * t188;
t196 = t186 * pkin(3) + t200;
t184 = qJ(1) + pkin(9);
t179 = sin(t184);
t180 = cos(t184);
t195 = g(1) * t180 + g(2) * t179;
t187 = sin(qJ(1));
t190 = cos(qJ(1));
t194 = -g(1) * t190 - g(2) * t187;
t193 = t190 * pkin(1) + t179 * pkin(7) + t203 * t180;
t192 = t187 * pkin(1) - t180 * pkin(7) + t203 * t179;
t167 = t179 * t185 - t180 * t197;
t169 = t179 * t197 + t180 * t185;
t191 = g(1) * t167 - g(2) * t169 + t188 * t201;
t170 = t179 * t198 - t180 * t188;
t168 = t179 * t188 + t180 * t198;
t166 = g(3) * t186 + t195 * t189;
t165 = t195 * t186 - t201;
t164 = -g(1) * t168 - g(2) * t170 + t185 * t201;
t1 = [0, t194, g(1) * t187 - g(2) * t190, t194 * pkin(1) - g(3) * t200, 0, 0, 0, 0, 0, -t166, t165, -g(1) * t179 + g(2) * t180, t166, -t165, -g(1) * t193 - g(2) * t192 - g(3) * (-t189 * qJ(4) + t196) 0, 0, 0, 0, 0, t164, t191, t164, -t166, -t191, -g(1) * (t179 * pkin(4) + t168 * pkin(5) + t167 * qJ(6) + t193) - g(2) * (-t180 * pkin(4) + t170 * pkin(5) - t169 * qJ(6) + t192) - g(3) * (t186 * pkin(8) + t196) + (-g(3) * (-pkin(5) * t185 + qJ(6) * t188 - qJ(4)) - t195 * pkin(8)) * t189;];
U_reg  = t1;
