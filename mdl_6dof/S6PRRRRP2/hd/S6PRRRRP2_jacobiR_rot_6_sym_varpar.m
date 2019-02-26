% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:43
% EndTime: 2019-02-26 20:15:43
% DurationCPUTime: 0.11s
% Computational Cost: add. (123->25), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
t186 = qJ(3) + qJ(4);
t185 = cos(t186);
t191 = sin(qJ(5));
t203 = t185 * t191;
t193 = cos(qJ(5));
t202 = t185 * t193;
t187 = sin(pkin(11));
t188 = sin(pkin(6));
t201 = t187 * t188;
t189 = cos(pkin(11));
t200 = t188 * t189;
t192 = sin(qJ(2));
t199 = t188 * t192;
t190 = cos(pkin(6));
t198 = t190 * t192;
t194 = cos(qJ(2));
t197 = t190 * t194;
t196 = t191 * t194;
t195 = t193 * t194;
t184 = sin(t186);
t182 = -t187 * t198 + t189 * t194;
t181 = t187 * t197 + t189 * t192;
t180 = t187 * t194 + t189 * t198;
t179 = t187 * t192 - t189 * t197;
t178 = t190 * t184 + t185 * t199;
t177 = -t184 * t199 + t190 * t185;
t176 = t177 * t193;
t175 = t177 * t191;
t174 = t182 * t185 + t184 * t201;
t173 = -t182 * t184 + t185 * t201;
t172 = t180 * t185 - t184 * t200;
t171 = -t180 * t184 - t185 * t200;
t170 = t173 * t193;
t169 = t173 * t191;
t168 = t171 * t193;
t167 = t171 * t191;
t1 = [0, -t181 * t202 + t182 * t191, t170, t170, -t174 * t191 + t181 * t193, 0; 0, -t179 * t202 + t180 * t191, t168, t168, -t172 * t191 + t179 * t193, 0; 0 (t185 * t195 + t191 * t192) * t188, t176, t176, -t178 * t191 - t188 * t195, 0; 0, -t181 * t184, t174, t174, 0, 0; 0, -t179 * t184, t172, t172, 0, 0; 0, t188 * t194 * t184, t178, t178, 0, 0; 0, -t181 * t203 - t182 * t193, t169, t169, t174 * t193 + t181 * t191, 0; 0, -t179 * t203 - t180 * t193, t167, t167, t172 * t193 + t179 * t191, 0; 0 (t185 * t196 - t192 * t193) * t188, t175, t175, t178 * t193 - t188 * t196, 0;];
JR_rot  = t1;
