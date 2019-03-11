% Calculate minimal parameter regressor of potential energy for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:34
% EndTime: 2019-03-09 03:49:35
% DurationCPUTime: 0.09s
% Computational Cost: add. (101->39), mult. (127->57), div. (0->0), fcn. (137->10), ass. (0->23)
t216 = sin(qJ(1));
t219 = cos(qJ(1));
t222 = g(1) * t219 + g(2) * t216;
t210 = pkin(10) + qJ(3);
t207 = sin(t210);
t208 = cos(t210);
t215 = sin(qJ(5));
t218 = cos(qJ(5));
t204 = t207 * t218 - t208 * t215;
t223 = g(3) * t204;
t221 = t207 * t215 + t208 * t218;
t212 = cos(pkin(10));
t220 = t212 * pkin(2) + pkin(3) * t208 + qJ(4) * t207 + pkin(1);
t217 = cos(qJ(6));
t214 = sin(qJ(6));
t213 = -pkin(7) - qJ(2);
t211 = sin(pkin(10));
t205 = g(1) * t216 - g(2) * t219;
t203 = t221 * t219;
t202 = t221 * t216;
t201 = -g(3) * t207 - t222 * t208;
t200 = -g(3) * t208 + t222 * t207;
t1 = [0, -t222, t205, -g(3) * t211 - t222 * t212, -g(3) * t212 + t222 * t211, -t205, -g(1) * (t219 * pkin(1) + t216 * qJ(2)) - g(2) * (t216 * pkin(1) - t219 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, t201, t200, t201, -t205, -t200, -g(3) * (t211 * pkin(2) + t207 * pkin(3) - t208 * qJ(4) + pkin(6)) + (-g(1) * t220 - g(2) * t213) * t219 + (g(1) * t213 - g(2) * t220) * t216, 0, 0, 0, 0, 0, -g(1) * t203 - g(2) * t202 - t223, g(3) * t221 - t222 * t204, 0, 0, 0, 0, 0, -g(1) * (t203 * t217 - t216 * t214) - g(2) * (t202 * t217 + t219 * t214) - t217 * t223, -g(1) * (-t203 * t214 - t216 * t217) - g(2) * (-t202 * t214 + t219 * t217) + t214 * t223;];
U_reg  = t1;
