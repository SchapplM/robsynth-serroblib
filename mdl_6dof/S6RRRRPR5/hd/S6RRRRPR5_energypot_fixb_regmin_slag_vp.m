% Calculate minimal parameter regressor of potential energy for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x35]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:16:16
% EndTime: 2019-03-09 22:16:16
% DurationCPUTime: 0.14s
% Computational Cost: add. (116->45), mult. (155->64), div. (0->0), fcn. (178->10), ass. (0->27)
t219 = qJ(2) + qJ(3);
t218 = cos(t219);
t226 = cos(qJ(2));
t237 = t226 * pkin(2) + pkin(3) * t218 + pkin(1);
t217 = sin(t219);
t235 = g(3) * t217;
t221 = sin(qJ(4));
t223 = sin(qJ(1));
t234 = t223 * t221;
t225 = cos(qJ(4));
t233 = t223 * t225;
t227 = cos(qJ(1));
t232 = t227 * t221;
t231 = t227 * t225;
t230 = g(1) * t227 + g(2) * t223;
t211 = t218 * t234 + t231;
t213 = t218 * t232 - t233;
t229 = g(1) * t213 + g(2) * t211 + t221 * t235;
t228 = -pkin(8) - pkin(7);
t224 = cos(qJ(6));
t222 = sin(qJ(2));
t220 = sin(qJ(6));
t214 = t218 * t231 + t234;
t212 = t218 * t233 - t232;
t210 = -g(3) * t218 + t230 * t217;
t209 = -g(1) * t214 - g(2) * t212 - t225 * t235;
t1 = [0, -t230, g(1) * t223 - g(2) * t227, 0, 0, 0, 0, 0, -g(3) * t222 - t230 * t226, -g(3) * t226 + t230 * t222, 0, 0, 0, 0, 0, -t230 * t218 - t235, t210, 0, 0, 0, 0, 0, t209, t229, t209, -t210, -t229, -g(1) * (t214 * pkin(4) + t213 * qJ(5) - t223 * t228 + t237 * t227) - g(2) * (t212 * pkin(4) + t211 * qJ(5) + t237 * t223 + t227 * t228) - g(3) * (t222 * pkin(2) - t218 * pkin(9) + pkin(6)) + (-g(3) * (pkin(4) * t225 + qJ(5) * t221 + pkin(3)) - t230 * pkin(9)) * t217, 0, 0, 0, 0, 0, -g(1) * (t213 * t220 + t214 * t224) - g(2) * (t211 * t220 + t212 * t224) - (t220 * t221 + t224 * t225) * t235, -g(1) * (t213 * t224 - t214 * t220) - g(2) * (t211 * t224 - t212 * t220) - (-t220 * t225 + t221 * t224) * t235;];
U_reg  = t1;
