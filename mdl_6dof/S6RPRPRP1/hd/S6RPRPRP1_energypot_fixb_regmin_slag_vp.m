% Calculate minimal parameter regressor of potential energy for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:02
% EndTime: 2019-03-09 03:03:02
% DurationCPUTime: 0.10s
% Computational Cost: add. (103->40), mult. (87->56), div. (0->0), fcn. (82->10), ass. (0->29)
t196 = cos(qJ(5));
t180 = t196 * pkin(5) + pkin(4);
t189 = qJ(3) + pkin(10);
t182 = sin(t189);
t184 = cos(t189);
t191 = -qJ(6) - pkin(8);
t212 = t180 * t184 - t182 * t191;
t211 = g(3) * t182;
t210 = qJ(2) + pkin(6);
t190 = qJ(1) + pkin(9);
t183 = sin(t190);
t193 = sin(qJ(5));
t207 = t183 * t193;
t206 = t183 * t196;
t185 = cos(t190);
t205 = t185 * t193;
t204 = t185 * t196;
t197 = cos(qJ(3));
t181 = t197 * pkin(3) + pkin(2);
t198 = cos(qJ(1));
t203 = t198 * pkin(1) + t185 * t181;
t194 = sin(qJ(3));
t202 = t194 * pkin(3) + t210;
t192 = -qJ(4) - pkin(7);
t195 = sin(qJ(1));
t201 = t195 * pkin(1) + t183 * t181 + t185 * t192;
t200 = g(1) * t185 + g(2) * t183;
t199 = -g(1) * t198 - g(2) * t195;
t1 = [0, t199, g(1) * t195 - g(2) * t198, t199 * pkin(1) - g(3) * t210, 0, 0, 0, 0, 0, -g(3) * t194 - t200 * t197, -g(3) * t197 + t200 * t194, -g(1) * t183 + g(2) * t185, -g(1) * (-t183 * t192 + t203) - g(2) * t201 - g(3) * t202, 0, 0, 0, 0, 0, -g(1) * (t184 * t204 + t207) - g(2) * (t184 * t206 - t205) - t196 * t211, -g(1) * (-t184 * t205 + t206) - g(2) * (-t184 * t207 - t204) + t193 * t211, g(3) * t184 - t200 * t182, -g(1) * (t212 * t185 + t203) - g(2) * (-pkin(5) * t205 + t201) - g(3) * (t182 * t180 + t184 * t191 + t202) + (-g(1) * (pkin(5) * t193 - t192) - g(2) * t212) * t183;];
U_reg  = t1;
