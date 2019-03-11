% Calculate minimal parameter regressor of potential energy for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:09:49
% EndTime: 2019-03-09 07:09:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (87->33), mult. (83->51), div. (0->0), fcn. (88->12), ass. (0->27)
t218 = pkin(11) + qJ(3);
t215 = qJ(4) + t218;
t211 = sin(t215);
t235 = g(3) * t211;
t219 = qJ(5) + qJ(6);
t216 = sin(t219);
t223 = sin(qJ(1));
t234 = t223 * t216;
t217 = cos(t219);
t233 = t223 * t217;
t222 = sin(qJ(5));
t232 = t223 * t222;
t224 = cos(qJ(5));
t231 = t223 * t224;
t225 = cos(qJ(1));
t230 = t225 * t216;
t229 = t225 * t217;
t228 = t225 * t222;
t227 = t225 * t224;
t226 = g(1) * t225 + g(2) * t223;
t221 = cos(pkin(11));
t220 = sin(pkin(11));
t214 = cos(t218);
t213 = sin(t218);
t212 = cos(t215);
t210 = g(1) * t223 - g(2) * t225;
t1 = [0, -t226, t210, -g(3) * t220 - t226 * t221, -g(3) * t221 + t226 * t220, -t210, -g(1) * (t225 * pkin(1) + t223 * qJ(2)) - g(2) * (t223 * pkin(1) - t225 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t213 - t226 * t214, -g(3) * t214 + t226 * t213, 0, 0, 0, 0, 0, -t226 * t212 - t235, -g(3) * t212 + t226 * t211, 0, 0, 0, 0, 0, -g(1) * (t212 * t227 + t232) - g(2) * (t212 * t231 - t228) - t224 * t235, -g(1) * (-t212 * t228 + t231) - g(2) * (-t212 * t232 - t227) + t222 * t235, 0, 0, 0, 0, 0, -g(1) * (t212 * t229 + t234) - g(2) * (t212 * t233 - t230) - t217 * t235, -g(1) * (-t212 * t230 + t233) - g(2) * (-t212 * t234 - t229) + t216 * t235;];
U_reg  = t1;
