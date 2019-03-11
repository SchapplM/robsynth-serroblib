% Calculate minimal parameter regressor of potential energy for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:42:38
% EndTime: 2019-03-09 02:42:38
% DurationCPUTime: 0.07s
% Computational Cost: add. (103->36), mult. (87->52), div. (0->0), fcn. (82->10), ass. (0->28)
t190 = qJ(3) + pkin(10);
t185 = cos(t190);
t210 = g(3) * t185;
t209 = qJ(2) + pkin(6);
t191 = qJ(1) + pkin(9);
t184 = sin(t191);
t193 = sin(qJ(6));
t208 = t184 * t193;
t196 = cos(qJ(6));
t207 = t184 * t196;
t186 = cos(t191);
t206 = t186 * t193;
t205 = t186 * t196;
t194 = sin(qJ(3));
t204 = t194 * pkin(3) + t209;
t197 = cos(qJ(3));
t182 = t197 * pkin(3) + pkin(2);
t192 = -qJ(4) - pkin(7);
t195 = sin(qJ(1));
t203 = t195 * pkin(1) + t184 * t182 + t186 * t192;
t202 = g(1) * t186 + g(2) * t184;
t198 = cos(qJ(1));
t201 = -g(1) * t198 - g(2) * t195;
t200 = t198 * pkin(1) + t186 * t182 - t184 * t192;
t183 = sin(t190);
t199 = pkin(4) * t185 + qJ(5) * t183;
t178 = -g(1) * t184 + g(2) * t186;
t1 = [0, t201, g(1) * t195 - g(2) * t198, t201 * pkin(1) - g(3) * t209, 0, 0, 0, 0, 0, -g(3) * t194 - t202 * t197, -g(3) * t197 + t202 * t194, t178, -g(1) * t200 - g(2) * t203 - g(3) * t204, t178, g(3) * t183 + t202 * t185, -t202 * t183 + t210, -g(1) * (t199 * t186 + t200) - g(2) * (t199 * t184 + t203) - g(3) * (t183 * pkin(4) - t185 * qJ(5) + t204) 0, 0, 0, 0, 0, -g(1) * (t183 * t206 + t207) - g(2) * (t183 * t208 - t205) + t193 * t210, -g(1) * (t183 * t205 - t208) - g(2) * (t183 * t207 + t206) + t196 * t210;];
U_reg  = t1;
