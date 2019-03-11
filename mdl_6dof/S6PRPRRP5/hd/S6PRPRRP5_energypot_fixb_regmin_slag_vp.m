% Calculate minimal parameter regressor of potential energy for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:42
% EndTime: 2019-03-08 20:16:42
% DurationCPUTime: 0.13s
% Computational Cost: add. (127->64), mult. (284->97), div. (0->0), fcn. (342->10), ass. (0->39)
t199 = sin(pkin(6));
t224 = pkin(7) * t199;
t198 = sin(pkin(10));
t200 = cos(pkin(10));
t208 = cos(qJ(2));
t201 = cos(pkin(6));
t205 = sin(qJ(2));
t217 = t201 * t205;
t185 = t198 * t208 + t200 * t217;
t203 = sin(qJ(5));
t223 = t185 * t203;
t187 = -t198 * t217 + t200 * t208;
t222 = t187 * t203;
t204 = sin(qJ(4));
t221 = t199 * t204;
t220 = t199 * t205;
t207 = cos(qJ(4));
t219 = t199 * t207;
t218 = t199 * t208;
t216 = t201 * t208;
t215 = t203 * t205;
t214 = pkin(2) * t220 + t201 * pkin(7) + qJ(1);
t184 = t198 * t205 - t200 * t216;
t213 = t198 * pkin(1) + t185 * pkin(2) + t184 * qJ(3);
t186 = t198 * t216 + t200 * t205;
t212 = t200 * pkin(1) + t187 * pkin(2) + t186 * qJ(3) + t198 * t224;
t177 = -t186 * t207 + t198 * t221;
t179 = t184 * t207 + t200 * t221;
t188 = t201 * t204 + t207 * t218;
t211 = g(1) * t177 - g(2) * t179 + g(3) * t188;
t210 = -g(1) * t186 - g(2) * t184 + g(3) * t218;
t209 = g(1) * t187 + g(2) * t185 + g(3) * t220;
t206 = cos(qJ(5));
t202 = -qJ(6) - pkin(9);
t194 = t206 * pkin(5) + pkin(4);
t189 = t201 * t207 - t204 * t218;
t180 = t184 * t204 - t200 * t219;
t178 = t186 * t204 + t198 * t219;
t1 = [-g(3) * qJ(1), 0, -t209, -t210, t209, t210, -g(1) * t212 - g(2) * (-t200 * t224 + t213) - g(3) * (-qJ(3) * t218 + t214) 0, 0, 0, 0, 0, -g(1) * t178 - g(2) * t180 - g(3) * t189, t211, 0, 0, 0, 0, 0, -g(1) * (t178 * t206 + t222) - g(2) * (t180 * t206 + t223) - g(3) * (t189 * t206 + t199 * t215) -g(1) * (-t178 * t203 + t187 * t206) - g(2) * (-t180 * t203 + t185 * t206) - g(3) * (-t189 * t203 + t206 * t220) -t211, -g(1) * (pkin(5) * t222 + t187 * pkin(8) - t177 * t202 + t178 * t194 + t212) - g(2) * (pkin(5) * t223 + t185 * pkin(8) + t179 * t202 + t180 * t194 + t213) - g(3) * (t201 * pkin(3) - t188 * t202 + t189 * t194 + t214) + (-g(1) * pkin(3) * t198 - g(3) * (pkin(5) * t215 + pkin(8) * t205 - qJ(3) * t208) - g(2) * (-pkin(3) - pkin(7)) * t200) * t199;];
U_reg  = t1;
