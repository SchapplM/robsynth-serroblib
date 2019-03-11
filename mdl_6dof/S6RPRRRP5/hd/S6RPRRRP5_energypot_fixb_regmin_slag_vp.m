% Calculate minimal parameter regressor of potential energy for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x32]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:12:24
% EndTime: 2019-03-09 06:12:24
% DurationCPUTime: 0.09s
% Computational Cost: add. (142->45), mult. (131->59), div. (0->0), fcn. (137->10), ass. (0->29)
t211 = pkin(10) + qJ(3);
t209 = qJ(4) + t211;
t205 = cos(t209);
t208 = cos(t211);
t213 = cos(pkin(10));
t226 = t213 * pkin(2) + pkin(3) * t208 + pkin(4) * t205 + pkin(1);
t204 = sin(t209);
t224 = g(3) * t204;
t214 = sin(qJ(5));
t215 = sin(qJ(1));
t223 = t215 * t214;
t216 = cos(qJ(5));
t222 = t215 * t216;
t217 = cos(qJ(1));
t221 = t217 * t214;
t220 = t217 * t216;
t219 = g(1) * t217 + g(2) * t215;
t198 = t205 * t223 + t220;
t200 = t205 * t221 - t222;
t218 = g(1) * t200 + g(2) * t198 + t214 * t224;
t212 = sin(pkin(10));
t210 = -pkin(8) - pkin(7) - qJ(2);
t207 = sin(t211);
t203 = g(1) * t215 - g(2) * t217;
t201 = t205 * t220 + t223;
t199 = t205 * t222 - t221;
t197 = -g(3) * t205 + t219 * t204;
t196 = -g(1) * t201 - g(2) * t199 - t216 * t224;
t1 = [0, -t219, t203, -g(3) * t212 - t219 * t213, -g(3) * t213 + t219 * t212, -t203, -g(1) * (t217 * pkin(1) + t215 * qJ(2)) - g(2) * (t215 * pkin(1) - t217 * qJ(2)) - g(3) * pkin(6), 0, 0, 0, 0, 0, -g(3) * t207 - t219 * t208, -g(3) * t208 + t219 * t207, 0, 0, 0, 0, 0, -t219 * t205 - t224, t197, 0, 0, 0, 0, 0, t196, t218, t196, -t197, -t218, -g(1) * (t201 * pkin(5) + t200 * qJ(6) - t215 * t210 + t226 * t217) - g(2) * (t199 * pkin(5) + t198 * qJ(6) + t217 * t210 + t226 * t215) - g(3) * (t212 * pkin(2) + pkin(3) * t207 - t205 * pkin(9) + pkin(6)) + (-g(3) * (pkin(5) * t216 + qJ(6) * t214 + pkin(4)) - t219 * pkin(9)) * t204;];
U_reg  = t1;
