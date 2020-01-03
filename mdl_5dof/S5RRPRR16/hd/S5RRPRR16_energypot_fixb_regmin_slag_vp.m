% Calculate minimal parameter regressor of potential energy for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR16_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:25
% EndTime: 2019-12-31 20:47:25
% DurationCPUTime: 0.08s
% Computational Cost: add. (73->44), mult. (177->75), div. (0->0), fcn. (216->10), ass. (0->29)
t188 = sin(pkin(5));
t192 = sin(qJ(2));
t208 = t188 * t192;
t193 = sin(qJ(1));
t207 = t188 * t193;
t196 = cos(qJ(2));
t206 = t188 * t196;
t197 = cos(qJ(1));
t205 = t188 * t197;
t204 = t193 * t192;
t203 = t193 * t196;
t202 = t197 * t192;
t201 = t197 * t196;
t200 = g(1) * t193 - g(2) * t197;
t189 = cos(pkin(5));
t183 = -t189 * t201 + t204;
t185 = t189 * t203 + t202;
t199 = -g(1) * t185 - g(2) * t183 + g(3) * t206;
t184 = t189 * t202 + t203;
t186 = -t189 * t204 + t201;
t198 = g(1) * t186 + g(2) * t184 + g(3) * t208;
t195 = cos(qJ(4));
t194 = cos(qJ(5));
t191 = sin(qJ(4));
t190 = sin(qJ(5));
t182 = t189 * t195 - t191 * t206;
t181 = t183 * t191 - t195 * t205;
t180 = t185 * t191 + t195 * t207;
t1 = [0, -g(1) * t197 - g(2) * t193, t200, 0, 0, 0, 0, 0, -t198, -t199, -g(3) * t189 - t200 * t188, t198, t199, -g(1) * (t197 * pkin(1) + t186 * pkin(2) + pkin(7) * t207 + t185 * qJ(3)) - g(2) * (t193 * pkin(1) + t184 * pkin(2) - pkin(7) * t205 + t183 * qJ(3)) - g(3) * (t189 * pkin(7) + pkin(6) + (pkin(2) * t192 - qJ(3) * t196) * t188), 0, 0, 0, 0, 0, -g(1) * t180 - g(2) * t181 - g(3) * t182, -g(1) * (t185 * t195 - t191 * t207) - g(2) * (t183 * t195 + t191 * t205) - g(3) * (-t189 * t191 - t195 * t206), 0, 0, 0, 0, 0, -g(1) * (t180 * t194 + t186 * t190) - g(2) * (t181 * t194 + t184 * t190) - g(3) * (t182 * t194 + t190 * t208), -g(1) * (-t180 * t190 + t186 * t194) - g(2) * (-t181 * t190 + t184 * t194) - g(3) * (-t182 * t190 + t194 * t208);];
U_reg = t1;
