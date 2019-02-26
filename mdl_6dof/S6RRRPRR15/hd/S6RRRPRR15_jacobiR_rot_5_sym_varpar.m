% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR15_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:25
% EndTime: 2019-02-26 22:24:25
% DurationCPUTime: 0.16s
% Computational Cost: add. (133->42), mult. (408->89), div. (0->0), fcn. (571->12), ass. (0->50)
t182 = sin(qJ(2));
t183 = sin(qJ(1));
t186 = cos(qJ(2));
t187 = cos(qJ(1));
t207 = cos(pkin(6));
t191 = t187 * t207;
t172 = t182 * t191 + t183 * t186;
t181 = sin(qJ(3));
t185 = cos(qJ(3));
t171 = t183 * t182 - t186 * t191;
t177 = sin(pkin(7));
t179 = cos(pkin(7));
t178 = sin(pkin(6));
t200 = t178 * t187;
t189 = t171 * t179 + t177 * t200;
t157 = t172 * t181 + t189 * t185;
t166 = -t171 * t177 + t179 * t200;
t180 = sin(qJ(5));
t184 = cos(qJ(5));
t212 = -t157 * t180 + t166 * t184;
t211 = t157 * t184 + t166 * t180;
t208 = -t172 * t185 + t189 * t181;
t204 = t177 * t178;
t203 = t177 * t180;
t202 = t177 * t184;
t201 = t178 * t183;
t199 = t179 * t181;
t198 = t179 * t185;
t197 = t181 * t182;
t196 = t181 * t186;
t195 = t182 * t185;
t194 = t185 * t186;
t193 = t182 * t204;
t192 = t183 * t207;
t190 = t207 * t177;
t173 = -t187 * t182 - t186 * t192;
t188 = t173 * t179 + t177 * t201;
t174 = -t182 * t192 + t187 * t186;
t170 = t207 * t179 - t186 * t204;
t169 = (t179 * t195 + t196) * t178;
t168 = -t173 * t177 + t179 * t201;
t165 = t181 * t190 + (t179 * t196 + t195) * t178;
t164 = -t185 * t190 + (-t179 * t194 + t197) * t178;
t163 = t173 * t181 + t174 * t198;
t162 = -t171 * t181 + t172 * t198;
t161 = t174 * t185 + t188 * t181;
t160 = t174 * t181 - t188 * t185;
t156 = t160 * t180 + t168 * t184;
t155 = t160 * t184 - t168 * t180;
t1 = [t212, t163 * t180 + t174 * t202, t161 * t180, 0, t155, 0; t156, t162 * t180 + t172 * t202, -t208 * t180, 0, t211, 0; 0, t169 * t180 + t184 * t193, t165 * t180, 0, t164 * t184 - t170 * t180, 0; -t211, t163 * t184 - t174 * t203, t161 * t184, 0, -t156, 0; t155, t162 * t184 - t172 * t203, -t208 * t184, 0, t212, 0; 0, t169 * t184 - t180 * t193, t165 * t184, 0, -t164 * t180 - t170 * t184, 0; t208, t173 * t185 - t174 * t199, -t160, 0, 0, 0; t161, -t171 * t185 - t172 * t199, -t157, 0, 0, 0; 0 (-t179 * t197 + t194) * t178, -t164, 0, 0, 0;];
JR_rot  = t1;
