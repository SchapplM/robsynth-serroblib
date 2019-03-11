% Calculate minimal parameter regressor of potential energy for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:31:56
% EndTime: 2019-03-09 08:31:56
% DurationCPUTime: 0.09s
% Computational Cost: add. (107->49), mult. (120->61), div. (0->0), fcn. (115->8), ass. (0->33)
t208 = cos(qJ(1));
t200 = qJ(2) + pkin(9);
t197 = cos(t200);
t219 = t197 * t208;
t196 = sin(t200);
t220 = qJ(4) * t196;
t224 = pkin(3) * t219 + t208 * t220;
t203 = sin(qJ(5));
t223 = pkin(5) * t203;
t222 = g(3) * t197;
t204 = sin(qJ(2));
t221 = t204 * pkin(2) + pkin(6);
t205 = sin(qJ(1));
t218 = t205 * t203;
t206 = cos(qJ(5));
t217 = t205 * t206;
t216 = t208 * t203;
t215 = t208 * t206;
t207 = cos(qJ(2));
t195 = t207 * pkin(2) + pkin(1);
t202 = -pkin(7) - qJ(3);
t214 = t205 * t195 + t208 * t202;
t213 = t196 * pkin(3) + t221;
t212 = t196 * t216;
t191 = t208 * t195;
t211 = -t205 * t202 + t191;
t210 = t214 + (pkin(3) * t197 + t220) * t205;
t209 = g(1) * t208 + g(2) * t205;
t201 = -qJ(6) - pkin(8);
t194 = t206 * pkin(5) + pkin(4);
t187 = g(1) * t205 - g(2) * t208;
t184 = g(3) * t196 + t209 * t197;
t1 = [0, -t209, t187, 0, 0, 0, 0, 0, -g(3) * t204 - t209 * t207, -g(3) * t207 + t209 * t204, -t187, -g(1) * t211 - g(2) * t214 - g(3) * t221, -t187, t184, -t209 * t196 + t222, -g(1) * (t211 + t224) - g(2) * t210 - g(3) * (-t197 * qJ(4) + t213) 0, 0, 0, 0, 0, -g(1) * (t212 + t217) - g(2) * (t196 * t218 - t215) + t203 * t222, -g(1) * (t196 * t215 - t218) - g(2) * (t196 * t217 + t216) + t206 * t222, -t184, -g(1) * (pkin(5) * t212 - t201 * t219 + t191 + t224) - g(2) * (-t208 * t194 + t210) - g(3) * (-t196 * t201 + (-qJ(4) - t223) * t197 + t213) + (-g(1) * (t194 - t202) - g(2) * (t196 * t223 - t197 * t201)) * t205;];
U_reg  = t1;
