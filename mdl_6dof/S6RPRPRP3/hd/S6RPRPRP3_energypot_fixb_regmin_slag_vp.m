% Calculate minimal parameter regressor of potential energy for
% S6RPRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:09:36
% EndTime: 2019-03-09 03:09:36
% DurationCPUTime: 0.10s
% Computational Cost: add. (169->58), mult. (154->81), div. (0->0), fcn. (161->10), ass. (0->35)
t199 = qJ(1) + pkin(9);
t192 = sin(t199);
t220 = g(2) * t192;
t203 = sin(qJ(3));
t219 = g(3) * t203;
t218 = qJ(2) + pkin(6);
t200 = sin(pkin(10));
t217 = t192 * t200;
t205 = cos(qJ(3));
t216 = t192 * t205;
t194 = cos(t199);
t215 = t194 * t205;
t214 = t200 * t205;
t201 = cos(pkin(10));
t213 = t201 * t205;
t204 = sin(qJ(1));
t212 = t204 * pkin(1) + t192 * pkin(2);
t206 = cos(qJ(1));
t211 = t206 * pkin(1) + t194 * pkin(2) + t192 * pkin(7);
t210 = g(1) * t194 + t220;
t209 = -g(1) * t206 - g(2) * t204;
t208 = pkin(3) * t205 + qJ(4) * t203;
t198 = pkin(10) + qJ(5);
t191 = sin(t198);
t193 = cos(t198);
t181 = t191 * t216 + t194 * t193;
t183 = t191 * t215 - t192 * t193;
t207 = g(1) * t183 + g(2) * t181 + t191 * t219;
t202 = -pkin(8) - qJ(4);
t190 = t201 * pkin(4) + pkin(3);
t185 = -g(3) * t205 + t210 * t203;
t184 = t192 * t191 + t193 * t215;
t182 = -t194 * t191 + t193 * t216;
t180 = -g(1) * t184 - g(2) * t182 - t193 * t219;
t1 = [0, t209, g(1) * t204 - g(2) * t206, t209 * pkin(1) - g(3) * t218, 0, 0, 0, 0, 0, -t210 * t205 - t219, t185, -g(1) * (t194 * t213 + t217) - g(2) * (t192 * t213 - t194 * t200) - t201 * t219, -g(1) * (t192 * t201 - t194 * t214) - g(2) * (-t192 * t214 - t194 * t201) + t200 * t219, -t185, -g(1) * (t208 * t194 + t211) - g(2) * (-t194 * pkin(7) + t208 * t192 + t212) - g(3) * (t203 * pkin(3) - t205 * qJ(4) + t218) 0, 0, 0, 0, 0, t180, t207, t180, -t185, -t207, -g(1) * (pkin(4) * t217 + t184 * pkin(5) + t183 * qJ(6) + t211) - g(2) * (t182 * pkin(5) + t181 * qJ(6) + t190 * t216 + t212) - g(3) * (t205 * t202 + t218) + (t202 * t220 - g(3) * (pkin(5) * t193 + qJ(6) * t191 + t190)) * t203 + (-g(1) * (t190 * t205 - t202 * t203) - g(2) * (-pkin(4) * t200 - pkin(7))) * t194;];
U_reg  = t1;
