% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:59:45
% EndTime: 2019-03-09 13:59:45
% DurationCPUTime: 0.11s
% Computational Cost: add. (69->37), mult. (141->58), div. (0->0), fcn. (164->10), ass. (0->22)
t222 = sin(qJ(1));
t226 = cos(qJ(1));
t230 = g(1) * t226 + g(2) * t222;
t220 = sin(qJ(4));
t221 = sin(qJ(2));
t224 = cos(qJ(4));
t225 = cos(qJ(2));
t229 = t225 * t220 - t221 * t224;
t231 = g(3) * t229;
t228 = t221 * t220 + t225 * t224;
t227 = pkin(2) * t225 + qJ(3) * t221 + pkin(1);
t223 = cos(qJ(5));
t219 = sin(qJ(5));
t218 = qJ(5) + qJ(6);
t217 = cos(t218);
t216 = sin(t218);
t215 = g(1) * t222 - g(2) * t226;
t213 = t228 * t226;
t212 = t228 * t222;
t211 = -g(3) * t221 - t230 * t225;
t210 = -g(3) * t225 + t230 * t221;
t1 = [0, -t230, t215, 0, 0, 0, 0, 0, t211, t210, t211, -t215, -t210, -g(3) * (t221 * pkin(2) - t225 * qJ(3) + pkin(6)) + (g(2) * pkin(7) - g(1) * t227) * t226 + (-g(1) * pkin(7) - g(2) * t227) * t222, 0, 0, 0, 0, 0, -g(1) * t213 - g(2) * t212 + t231, g(3) * t228 + t230 * t229, 0, 0, 0, 0, 0, -g(1) * (t213 * t223 - t222 * t219) - g(2) * (t212 * t223 + t226 * t219) + t223 * t231, -g(1) * (-t213 * t219 - t222 * t223) - g(2) * (-t212 * t219 + t226 * t223) - t219 * t231, 0, 0, 0, 0, 0, -g(1) * (t213 * t217 - t222 * t216) - g(2) * (t212 * t217 + t226 * t216) + t217 * t231, -g(1) * (-t213 * t216 - t222 * t217) - g(2) * (-t212 * t216 + t226 * t217) - t216 * t231;];
U_reg  = t1;
