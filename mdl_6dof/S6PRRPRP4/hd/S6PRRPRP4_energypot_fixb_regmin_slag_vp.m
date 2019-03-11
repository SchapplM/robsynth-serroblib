% Calculate minimal parameter regressor of potential energy for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:43:58
% EndTime: 2019-03-08 21:43:58
% DurationCPUTime: 0.12s
% Computational Cost: add. (161->65), mult. (367->92), div. (0->0), fcn. (451->10), ass. (0->37)
t212 = sin(pkin(6));
t236 = pkin(7) * t212;
t219 = cos(qJ(5));
t235 = t219 * pkin(5) + pkin(4) + pkin(8);
t217 = sin(qJ(3));
t234 = t212 * t217;
t218 = sin(qJ(2));
t233 = t212 * t218;
t220 = cos(qJ(3));
t232 = t212 * t220;
t221 = cos(qJ(2));
t231 = t212 * t221;
t214 = cos(pkin(6));
t230 = t214 * t218;
t229 = t214 * t221;
t216 = sin(qJ(5));
t228 = pkin(5) * t216 + qJ(4);
t200 = t214 * t217 + t218 * t232;
t227 = pkin(2) * t233 + t200 * pkin(3) + t214 * pkin(7) + qJ(1);
t211 = sin(pkin(10));
t213 = cos(pkin(10));
t198 = -t211 * t230 + t213 * t221;
t191 = t198 * t220 + t211 * t234;
t226 = t213 * pkin(1) + t198 * pkin(2) + t191 * pkin(3) + t211 * t236;
t196 = t211 * t221 + t213 * t230;
t189 = t196 * t220 - t213 * t234;
t225 = t211 * pkin(1) + t196 * pkin(2) + t189 * pkin(3) - t213 * t236;
t188 = t196 * t217 + t213 * t232;
t190 = t198 * t217 - t211 * t232;
t199 = -t214 * t220 + t217 * t233;
t224 = g(1) * t190 + g(2) * t188 + g(3) * t199;
t223 = g(1) * t191 + g(2) * t189 + g(3) * t200;
t195 = t211 * t218 - t213 * t229;
t197 = t211 * t229 + t213 * t218;
t222 = -g(1) * t197 - g(2) * t195 + g(3) * t231;
t215 = -qJ(6) - pkin(9);
t1 = [-g(3) * qJ(1), 0, -g(1) * t198 - g(2) * t196 - g(3) * t233, -t222, 0, 0, 0, 0, 0, -t223, t224, t222, t223, -t224, -g(1) * (t197 * pkin(8) + t190 * qJ(4) + t226) - g(2) * (t195 * pkin(8) + t188 * qJ(4) + t225) - g(3) * (-pkin(8) * t231 + t199 * qJ(4) + t227) 0, 0, 0, 0, 0, -g(1) * (t190 * t216 + t197 * t219) - g(2) * (t188 * t216 + t195 * t219) - g(3) * (t199 * t216 - t219 * t231) -g(1) * (t190 * t219 - t197 * t216) - g(2) * (t188 * t219 - t195 * t216) - g(3) * (t199 * t219 + t216 * t231) -t223, -g(1) * (t228 * t190 - t191 * t215 + t235 * t197 + t226) - g(2) * (t228 * t188 - t189 * t215 + t235 * t195 + t225) - g(3) * (t228 * t199 - t200 * t215 - t235 * t231 + t227);];
U_reg  = t1;
