% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR6_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobiR_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:52
% EndTime: 2019-02-26 20:21:52
% DurationCPUTime: 0.17s
% Computational Cost: add. (177->51), mult. (541->118), div. (0->0), fcn. (739->14), ass. (0->55)
t181 = sin(pkin(8));
t182 = sin(pkin(7));
t217 = t181 * t182;
t183 = sin(pkin(6));
t216 = t182 * t183;
t185 = cos(pkin(8));
t215 = t182 * t185;
t187 = cos(pkin(6));
t214 = t182 * t187;
t186 = cos(pkin(7));
t213 = t183 * t186;
t188 = sin(qJ(4));
t212 = t185 * t188;
t191 = cos(qJ(4));
t211 = t185 * t191;
t189 = sin(qJ(3));
t210 = t186 * t189;
t192 = cos(qJ(3));
t209 = t186 * t192;
t190 = sin(qJ(2));
t208 = t187 * t190;
t193 = cos(qJ(2));
t207 = t187 * t193;
t206 = t189 * t190;
t205 = t189 * t193;
t204 = t190 * t192;
t203 = t192 * t193;
t184 = cos(pkin(14));
t202 = t184 * t216;
t201 = t190 * t216;
t180 = sin(pkin(14));
t174 = -t180 * t190 + t184 * t207;
t175 = t180 * t193 + t184 * t208;
t159 = -t175 * t189 + (t174 * t186 - t202) * t192;
t200 = t159 * t185 + (-t174 * t182 - t184 * t213) * t181;
t177 = -t180 * t208 + t184 * t193;
t176 = -t180 * t207 - t184 * t190;
t195 = t176 * t186 + t180 * t216;
t161 = -t177 * t189 + t195 * t192;
t199 = t161 * t185 + (-t176 * t182 + t180 * t213) * t181;
t167 = t192 * t214 + (t186 * t203 - t206) * t183;
t198 = t167 * t185 + (t187 * t186 - t193 * t216) * t181;
t163 = -t174 * t189 - t175 * t209;
t197 = t163 * t185 + t175 * t217;
t165 = -t176 * t189 - t177 * t209;
t196 = t165 * t185 + t177 * t217;
t171 = (-t186 * t204 - t205) * t183;
t194 = t171 * t185 + t181 * t201;
t172 = (-t186 * t206 + t203) * t183;
t168 = t189 * t214 + (t186 * t205 + t204) * t183;
t166 = t176 * t192 - t177 * t210;
t164 = t174 * t192 - t175 * t210;
t162 = t177 * t192 + t195 * t189;
t160 = t174 * t210 + t175 * t192 - t189 * t202;
t1 = [0, t166 * t191 + t196 * t188, t161 * t191 - t162 * t212, -t162 * t188 + t199 * t191, 0, 0; 0, t164 * t191 + t197 * t188, t159 * t191 - t160 * t212, -t160 * t188 + t200 * t191, 0, 0; 0, t172 * t191 + t194 * t188, t167 * t191 - t168 * t212, -t168 * t188 + t198 * t191, 0, 0; 0, -t166 * t188 + t196 * t191, -t161 * t188 - t162 * t211, -t162 * t191 - t199 * t188, 0, 0; 0, -t164 * t188 + t197 * t191, -t159 * t188 - t160 * t211, -t160 * t191 - t200 * t188, 0, 0; 0, -t172 * t188 + t194 * t191, -t167 * t188 - t168 * t211, -t168 * t191 - t198 * t188, 0, 0; 0, -t165 * t181 + t177 * t215, t162 * t181, 0, 0, 0; 0, -t163 * t181 + t175 * t215, t160 * t181, 0, 0, 0; 0, -t171 * t181 + t185 * t201, t168 * t181, 0, 0, 0;];
JR_rot  = t1;
