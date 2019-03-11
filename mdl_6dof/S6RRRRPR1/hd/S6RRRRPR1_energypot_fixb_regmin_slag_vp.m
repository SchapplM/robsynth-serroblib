% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:54:54
% EndTime: 2019-03-09 21:54:54
% DurationCPUTime: 0.06s
% Computational Cost: add. (87->32), mult. (72->44), div. (0->0), fcn. (73->12), ass. (0->24)
t220 = qJ(2) + qJ(3);
t219 = qJ(4) + t220;
t213 = pkin(11) + t219;
t232 = g(3) * sin(t213);
t221 = sin(qJ(6));
t223 = sin(qJ(1));
t231 = t223 * t221;
t224 = cos(qJ(6));
t230 = t223 * t224;
t226 = cos(qJ(1));
t229 = t226 * t221;
t228 = t226 * t224;
t227 = g(1) * t226 + g(2) * t223;
t225 = cos(qJ(2));
t222 = sin(qJ(2));
t218 = -qJ(5) - pkin(9) - pkin(8) - pkin(7);
t217 = cos(t220);
t216 = sin(t220);
t215 = cos(t219);
t214 = sin(t219);
t212 = cos(t213);
t210 = g(1) * t223 - g(2) * t226;
t209 = t225 * pkin(2) + pkin(3) * t217 + pkin(4) * t215 + pkin(1);
t1 = [0, -t227, t210, 0, 0, 0, 0, 0, -g(3) * t222 - t227 * t225, -g(3) * t225 + t227 * t222, 0, 0, 0, 0, 0, -g(3) * t216 - t227 * t217, -g(3) * t217 + t227 * t216, 0, 0, 0, 0, 0, -g(3) * t214 - t227 * t215, -g(3) * t215 + t227 * t214, -t210, -g(1) * (t226 * t209 - t223 * t218) - g(2) * (t223 * t209 + t226 * t218) - g(3) * (t222 * pkin(2) + pkin(3) * t216 + pkin(4) * t214 + pkin(6)) 0, 0, 0, 0, 0, -g(1) * (t212 * t228 + t231) - g(2) * (t212 * t230 - t229) - t224 * t232, -g(1) * (-t212 * t229 + t230) - g(2) * (-t212 * t231 - t228) + t221 * t232;];
U_reg  = t1;
