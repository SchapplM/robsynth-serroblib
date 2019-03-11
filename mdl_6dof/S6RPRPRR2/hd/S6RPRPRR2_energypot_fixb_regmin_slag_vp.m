% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:38:57
% EndTime: 2019-03-09 03:38:57
% DurationCPUTime: 0.07s
% Computational Cost: add. (86->36), mult. (73->54), div. (0->0), fcn. (75->12), ass. (0->29)
t206 = qJ(3) + pkin(11);
t227 = g(3) * sin(t206);
t226 = qJ(2) + pkin(6);
t207 = qJ(1) + pkin(10);
t201 = sin(t207);
t208 = qJ(5) + qJ(6);
t204 = sin(t208);
t225 = t201 * t204;
t205 = cos(t208);
t224 = t201 * t205;
t210 = sin(qJ(5));
t223 = t201 * t210;
t213 = cos(qJ(5));
t222 = t201 * t213;
t203 = cos(t207);
t221 = t203 * t204;
t220 = t203 * t205;
t219 = t203 * t210;
t218 = t203 * t213;
t217 = g(1) * t203 + g(2) * t201;
t212 = sin(qJ(1));
t215 = cos(qJ(1));
t216 = -g(1) * t215 - g(2) * t212;
t214 = cos(qJ(3));
t211 = sin(qJ(3));
t209 = -qJ(4) - pkin(7);
t202 = cos(t206);
t199 = t214 * pkin(3) + pkin(2);
t1 = [0, t216, g(1) * t212 - g(2) * t215, t216 * pkin(1) - g(3) * t226, 0, 0, 0, 0, 0, -g(3) * t211 - t217 * t214, -g(3) * t214 + t217 * t211, -g(1) * t201 + g(2) * t203, -g(1) * (t215 * pkin(1) + t203 * t199 - t201 * t209) - g(2) * (t212 * pkin(1) + t201 * t199 + t203 * t209) - g(3) * (t211 * pkin(3) + t226) 0, 0, 0, 0, 0, -g(1) * (t202 * t218 + t223) - g(2) * (t202 * t222 - t219) - t213 * t227, -g(1) * (-t202 * t219 + t222) - g(2) * (-t202 * t223 - t218) + t210 * t227, 0, 0, 0, 0, 0, -g(1) * (t202 * t220 + t225) - g(2) * (t202 * t224 - t221) - t205 * t227, -g(1) * (-t202 * t221 + t224) - g(2) * (-t202 * t225 - t220) + t204 * t227;];
U_reg  = t1;
