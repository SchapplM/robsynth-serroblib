% Calculate minimal parameter regressor of potential energy for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:52:28
% EndTime: 2019-03-09 09:52:28
% DurationCPUTime: 0.10s
% Computational Cost: add. (156->49), mult. (191->62), div. (0->0), fcn. (204->8), ass. (0->36)
t217 = qJ(2) + pkin(9);
t215 = cos(t217);
t241 = pkin(3) * t215;
t220 = sin(qJ(2));
t240 = t220 * pkin(2) + pkin(6);
t214 = sin(t217);
t219 = sin(qJ(4));
t239 = t214 * t219;
t221 = sin(qJ(1));
t238 = t214 * t221;
t222 = cos(qJ(4));
t237 = t214 * t222;
t224 = cos(qJ(1));
t236 = t214 * t224;
t235 = t221 * t219;
t234 = t221 * t222;
t233 = t224 * t219;
t232 = t224 * t222;
t223 = cos(qJ(2));
t213 = t223 * pkin(2) + pkin(1);
t218 = -pkin(7) - qJ(3);
t231 = t221 * t213 + t224 * t218;
t230 = t224 * t213 - t221 * t218;
t229 = g(1) * t224 + g(2) * t221;
t228 = t214 * pkin(3) + pkin(4) * t237 + qJ(5) * t239 + t240;
t197 = t215 * t235 + t232;
t198 = t215 * t234 - t233;
t227 = t198 * pkin(4) + pkin(8) * t238 + t197 * qJ(5) + t221 * t241 + t231;
t199 = t215 * t233 - t234;
t226 = g(1) * t199 + g(2) * t197 + g(3) * t239;
t200 = t215 * t232 + t235;
t225 = t200 * pkin(4) + pkin(8) * t236 + t199 * qJ(5) + t224 * t241 + t230;
t202 = g(1) * t221 - g(2) * t224;
t194 = -g(3) * t215 + t229 * t214;
t193 = -g(1) * t200 - g(2) * t198 - g(3) * t237;
t1 = [0, -t229, t202, 0, 0, 0, 0, 0, -g(3) * t220 - t229 * t223, -g(3) * t223 + t229 * t220, -t202, -g(1) * t230 - g(2) * t231 - g(3) * t240, 0, 0, 0, 0, 0, t193, t226, t193, -t194, -t226, -g(1) * t225 - g(2) * t227 - g(3) * (-t215 * pkin(8) + t228) t193, -t226, t194, -g(1) * (t200 * pkin(5) - qJ(6) * t236 + t225) - g(2) * (t198 * pkin(5) - qJ(6) * t238 + t227) - g(3) * (pkin(5) * t237 + (-pkin(8) + qJ(6)) * t215 + t228);];
U_reg  = t1;
