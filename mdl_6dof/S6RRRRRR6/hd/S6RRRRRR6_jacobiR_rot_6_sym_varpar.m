% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:04
% EndTime: 2019-02-26 22:50:05
% DurationCPUTime: 0.17s
% Computational Cost: add. (249->37), mult. (326->62), div. (0->0), fcn. (477->10), ass. (0->46)
t189 = cos(pkin(6));
t190 = sin(qJ(2));
t193 = cos(qJ(1));
t195 = t193 * t190;
t191 = sin(qJ(1));
t192 = cos(qJ(2));
t196 = t191 * t192;
t177 = t189 * t195 + t196;
t187 = qJ(3) + qJ(4);
t183 = sin(t187);
t185 = cos(t187);
t188 = sin(pkin(6));
t198 = t188 * t193;
t170 = -t177 * t185 + t183 * t198;
t194 = t193 * t192;
t197 = t191 * t190;
t176 = -t189 * t194 + t197;
t186 = qJ(5) + qJ(6);
t182 = sin(t186);
t184 = cos(t186);
t160 = t170 * t182 + t176 * t184;
t161 = t170 * t184 - t176 * t182;
t168 = -t177 * t183 - t185 * t198;
t209 = t168 * t182;
t179 = -t189 * t197 + t194;
t200 = t188 * t191;
t171 = t179 * t183 - t185 * t200;
t208 = t171 * t182;
t201 = t188 * t190;
t174 = -t183 * t201 + t189 * t185;
t207 = t174 * t182;
t204 = t182 * t185;
t203 = t184 * t185;
t202 = t185 * t192;
t199 = t188 * t192;
t178 = t189 * t196 + t195;
t175 = t189 * t183 + t185 * t201;
t173 = t174 * t184;
t172 = t179 * t185 + t183 * t200;
t167 = -t175 * t184 + t182 * t199;
t166 = -t175 * t182 - t184 * t199;
t165 = t171 * t184;
t164 = t168 * t184;
t163 = t172 * t184 + t178 * t182;
t162 = -t172 * t182 + t178 * t184;
t1 = [t161, -t178 * t203 + t179 * t182, -t165, -t165, t162, t162; t163, -t176 * t203 + t177 * t182, t164, t164, t160, t160; 0 (t182 * t190 + t184 * t202) * t188, t173, t173, t166, t166; -t160, t178 * t204 + t179 * t184, t208, t208, -t163, -t163; t162, t176 * t204 + t177 * t184, -t209, -t209, t161, t161; 0 (-t182 * t202 + t184 * t190) * t188, -t207, -t207, t167, t167; t168, -t178 * t183, t172, t172, 0, 0; t171, -t176 * t183, -t170, -t170, 0, 0; 0, t183 * t199, t175, t175, 0, 0;];
JR_rot  = t1;
