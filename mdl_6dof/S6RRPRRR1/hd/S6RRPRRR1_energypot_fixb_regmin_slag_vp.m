% Calculate minimal parameter regressor of potential energy for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x33]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:14:28
% EndTime: 2019-03-09 13:14:28
% DurationCPUTime: 0.05s
% Computational Cost: add. (80->27), mult. (66->39), div. (0->0), fcn. (67->10), ass. (0->22)
t218 = qJ(2) + pkin(11) + qJ(4);
t216 = qJ(5) + t218;
t212 = sin(t216);
t231 = g(3) * t212;
t220 = sin(qJ(6));
t222 = sin(qJ(1));
t230 = t222 * t220;
t223 = cos(qJ(6));
t229 = t222 * t223;
t225 = cos(qJ(1));
t228 = t225 * t220;
t227 = t225 * t223;
t226 = g(1) * t225 + g(2) * t222;
t224 = cos(qJ(2));
t221 = sin(qJ(2));
t219 = -pkin(7) - qJ(3);
t217 = t224 * pkin(2) + pkin(1);
t215 = cos(t218);
t214 = sin(t218);
t213 = cos(t216);
t211 = g(1) * t222 - g(2) * t225;
t1 = [0, -t226, t211, 0, 0, 0, 0, 0, -g(3) * t221 - t226 * t224, -g(3) * t224 + t226 * t221, -t211, -g(1) * (t225 * t217 - t222 * t219) - g(2) * (t222 * t217 + t225 * t219) - g(3) * (t221 * pkin(2) + pkin(6)) 0, 0, 0, 0, 0, -g(3) * t214 - t226 * t215, -g(3) * t215 + t226 * t214, 0, 0, 0, 0, 0, -t226 * t213 - t231, -g(3) * t213 + t226 * t212, 0, 0, 0, 0, 0, -g(1) * (t213 * t227 + t230) - g(2) * (t213 * t229 - t228) - t223 * t231, -g(1) * (-t213 * t228 + t229) - g(2) * (-t213 * t230 - t227) + t220 * t231;];
U_reg  = t1;
