% Calculate minimal parameter regressor of potential energy for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% U_reg [1x38]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR10V2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:51:47
% EndTime: 2019-04-11 14:51:47
% DurationCPUTime: 0.08s
% Computational Cost: add. (98->36), mult. (142->66), div. (0->0), fcn. (176->12), ass. (0->30)
t185 = qJ(2) + qJ(3);
t183 = sin(t185);
t188 = sin(qJ(4));
t204 = t183 * t188;
t190 = sin(qJ(1));
t203 = t183 * t190;
t193 = cos(qJ(4));
t202 = t183 * t193;
t195 = cos(qJ(1));
t201 = t183 * t195;
t200 = t190 * t188;
t199 = t190 * t193;
t198 = t195 * t188;
t197 = t195 * t193;
t196 = g(1) * t195 + g(2) * t190;
t194 = cos(qJ(2));
t192 = cos(qJ(5));
t191 = cos(qJ(6));
t189 = sin(qJ(2));
t187 = sin(qJ(5));
t186 = sin(qJ(6));
t184 = cos(t185);
t182 = t184 * t197 + t200;
t181 = t184 * t198 - t199;
t180 = t184 * t199 - t198;
t179 = t184 * t200 + t197;
t178 = -t184 * t187 + t192 * t202;
t177 = t182 * t192 + t187 * t201;
t176 = t180 * t192 + t187 * t203;
t1 = [0, -t196, g(1) * t190 - g(2) * t195, 0, 0, 0, 0, 0, -g(3) * t189 - t196 * t194, -g(3) * t194 + t196 * t189, 0, 0, 0, 0, 0, -g(3) * t183 - t196 * t184, -g(3) * t184 + t196 * t183, 0, 0, 0, 0, 0, -g(1) * t182 - g(2) * t180 - g(3) * t202, g(1) * t181 + g(2) * t179 + g(3) * t204, 0, 0, 0, 0, 0, -g(1) * t177 - g(2) * t176 - g(3) * t178, -g(1) * (-t182 * t187 + t192 * t201) - g(2) * (-t180 * t187 + t192 * t203) - g(3) * (-t184 * t192 - t187 * t202) 0, 0, 0, 0, 0, -g(1) * (t177 * t191 + t181 * t186) - g(2) * (t176 * t191 + t179 * t186) - g(3) * (t178 * t191 + t186 * t204) -g(1) * (-t177 * t186 + t181 * t191) - g(2) * (-t176 * t186 + t179 * t191) - g(3) * (-t178 * t186 + t191 * t204);];
U_reg  = t1;
