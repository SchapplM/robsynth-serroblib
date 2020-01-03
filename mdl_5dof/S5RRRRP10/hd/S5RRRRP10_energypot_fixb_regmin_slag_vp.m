% Calculate minimal parameter regressor of potential energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:20
% EndTime: 2019-12-31 22:12:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (97->52), mult. (222->84), div. (0->0), fcn. (274->10), ass. (0->35)
t202 = cos(pkin(5));
t210 = cos(qJ(2));
t211 = cos(qJ(1));
t214 = t211 * t210;
t206 = sin(qJ(2));
t207 = sin(qJ(1));
t217 = t207 * t206;
t193 = -t202 * t214 + t217;
t204 = sin(qJ(4));
t223 = t193 * t204;
t215 = t211 * t206;
t216 = t207 * t210;
t195 = t202 * t216 + t215;
t222 = t195 * t204;
t201 = sin(pkin(5));
t221 = t201 * t206;
t209 = cos(qJ(3));
t220 = t201 * t209;
t219 = t201 * t210;
t218 = t201 * t211;
t213 = g(1) * t207 - g(2) * t211;
t194 = t202 * t215 + t216;
t205 = sin(qJ(3));
t187 = t194 * t205 + t209 * t218;
t196 = -t202 * t217 + t214;
t189 = t196 * t205 - t207 * t220;
t191 = -t202 * t209 + t205 * t221;
t212 = g(1) * t189 + g(2) * t187 + g(3) * t191;
t208 = cos(qJ(4));
t203 = -qJ(5) - pkin(9);
t200 = t208 * pkin(4) + pkin(3);
t192 = t202 * t205 + t206 * t220;
t190 = t207 * t201 * t205 + t196 * t209;
t188 = t194 * t209 - t205 * t218;
t1 = [0, -g(1) * t211 - g(2) * t207, t213, 0, 0, 0, 0, 0, -g(1) * t196 - g(2) * t194 - g(3) * t221, g(1) * t195 + g(2) * t193 - g(3) * t219, 0, 0, 0, 0, 0, -g(1) * t190 - g(2) * t188 - g(3) * t192, t212, 0, 0, 0, 0, 0, -g(1) * (t190 * t208 + t222) - g(2) * (t188 * t208 + t223) - g(3) * (t192 * t208 - t204 * t219), -g(1) * (-t190 * t204 + t195 * t208) - g(2) * (-t188 * t204 + t193 * t208) - g(3) * (-t192 * t204 - t208 * t219), -t212, -g(1) * (t211 * pkin(1) + t196 * pkin(2) + pkin(4) * t222 + t195 * pkin(8) - t189 * t203 + t190 * t200) - g(2) * (t207 * pkin(1) + t194 * pkin(2) + pkin(4) * t223 + t193 * pkin(8) - t187 * t203 + t188 * t200) - g(3) * (t202 * pkin(7) - t191 * t203 + t192 * t200 + pkin(6)) + (-g(3) * (pkin(2) * t206 + (-pkin(4) * t204 - pkin(8)) * t210) - t213 * pkin(7)) * t201;];
U_reg = t1;
