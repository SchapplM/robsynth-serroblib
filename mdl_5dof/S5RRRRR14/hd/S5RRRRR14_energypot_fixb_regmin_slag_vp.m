% Calculate minimal parameter regressor of potential energy for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR14_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR14_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:31
% EndTime: 2024-09-27 18:43:31
% DurationCPUTime: 0.02s
% Computational Cost: add. (140->43), mult. (80->65), div. (0->0), fcn. (80->20), ass. (0->34)
t209 = g(3) * sin(pkin(5));
t202 = cos(pkin(5));
t203 = sin(qJ(3));
t208 = t202 * t203;
t205 = cos(qJ(3));
t207 = t202 * t205;
t199 = qJ(3) + qJ(4);
t193 = pkin(5) - t199;
t192 = pkin(5) + t199;
t206 = cos(qJ(1));
t204 = sin(qJ(1));
t200 = qJ(1) + qJ(2);
t198 = qJ(5) + t199;
t197 = cos(t200);
t196 = cos(t199);
t195 = sin(t200);
t194 = sin(t199);
t191 = cos(t198);
t190 = sin(t198);
t189 = -qJ(5) + t193;
t188 = qJ(5) + t192;
t187 = cos(t192);
t186 = sin(t193);
t185 = cos(t188);
t184 = sin(t189);
t183 = cos(t193) / 0.2e1;
t182 = sin(t192) / 0.2e1;
t181 = cos(t189) / 0.2e1;
t180 = sin(t188) / 0.2e1;
t179 = t183 + t187 / 0.2e1;
t178 = t182 - t186 / 0.2e1;
t177 = t181 + t185 / 0.2e1;
t176 = t180 - t184 / 0.2e1;
t1 = [0, -g(1) * t206 - g(2) * t204, g(1) * t204 - g(2) * t206, 0, -g(1) * t197 - g(2) * t195, g(1) * t195 - g(2) * t197, 0, 0, 0, 0, 0, -g(1) * (-t195 * t208 + t197 * t205) - g(2) * (t195 * t205 + t197 * t208) - t203 * t209, -g(1) * (-t195 * t207 - t197 * t203) - g(2) * (-t195 * t203 + t197 * t207) - t205 * t209, 0, 0, 0, 0, 0, -g(1) * (-t195 * t178 + t197 * t196) - g(2) * (t197 * t178 + t195 * t196) - g(3) * (t183 - t187 / 0.2e1), -g(1) * (-t195 * t179 - t197 * t194) - g(2) * (t197 * t179 - t195 * t194) - g(3) * (t182 + t186 / 0.2e1), 0, 0, 0, 0, 0, -g(1) * (-t195 * t176 + t197 * t191) - g(2) * (t197 * t176 + t195 * t191) - g(3) * (t181 - t185 / 0.2e1), -g(1) * (-t195 * t177 - t197 * t190) - g(2) * (t197 * t177 - t195 * t190) - g(3) * (t180 + t184 / 0.2e1);];
U_reg = t1;
