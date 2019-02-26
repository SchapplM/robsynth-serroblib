% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR12_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:59
% EndTime: 2019-02-26 22:36:59
% DurationCPUTime: 0.18s
% Computational Cost: add. (174->43), mult. (408->90), div. (0->0), fcn. (571->12), ass. (0->51)
t194 = sin(qJ(2));
t195 = sin(qJ(1));
t197 = cos(qJ(2));
t198 = cos(qJ(1));
t217 = cos(pkin(6));
t201 = t198 * t217;
t181 = t194 * t201 + t195 * t197;
t193 = sin(qJ(3));
t196 = cos(qJ(3));
t180 = t195 * t194 - t197 * t201;
t190 = sin(pkin(7));
t192 = cos(pkin(7));
t191 = sin(pkin(6));
t211 = t191 * t198;
t199 = t180 * t192 + t190 * t211;
t167 = -t181 * t196 + t199 * t193;
t174 = -t180 * t190 + t192 * t211;
t189 = qJ(4) + pkin(13);
t187 = sin(t189);
t188 = cos(t189);
t221 = t167 * t187 - t174 * t188;
t220 = t167 * t188 + t174 * t187;
t215 = t187 * t190;
t214 = t188 * t190;
t213 = t190 * t191;
t212 = t191 * t195;
t210 = t192 * t193;
t209 = t192 * t196;
t208 = t193 * t194;
t207 = t193 * t197;
t206 = t194 * t196;
t205 = t196 * t197;
t204 = t194 * t213;
t203 = t190 * t212;
t202 = t195 * t217;
t200 = t217 * t190;
t165 = -t181 * t193 - t199 * t196;
t183 = -t194 * t202 + t198 * t197;
t182 = -t198 * t194 - t197 * t202;
t179 = t217 * t192 - t197 * t213;
t178 = (-t192 * t208 + t205) * t191;
t176 = -t182 * t190 + t192 * t212;
t173 = t193 * t200 + (t192 * t207 + t206) * t191;
t172 = t196 * t200 + (t192 * t205 - t208) * t191;
t171 = t182 * t196 - t183 * t210;
t170 = -t180 * t196 - t181 * t210;
t169 = t183 * t196 + (t182 * t192 + t203) * t193;
t168 = -t182 * t209 + t183 * t193 - t196 * t203;
t164 = t169 * t188 + t176 * t187;
t163 = -t169 * t187 + t176 * t188;
t1 = [t220, t171 * t188 + t183 * t215, -t168 * t188, t163, 0, 0; t164, t170 * t188 + t181 * t215, t165 * t188, t221, 0, 0; 0, t178 * t188 + t187 * t204, t172 * t188, -t173 * t187 + t179 * t188, 0, 0; -t221, -t171 * t187 + t183 * t214, t168 * t187, -t164, 0, 0; t163, -t170 * t187 + t181 * t214, -t165 * t187, t220, 0, 0; 0, -t178 * t187 + t188 * t204, -t172 * t187, -t173 * t188 - t179 * t187, 0, 0; t165, t182 * t193 + t183 * t209, t169, 0, 0, 0; t168, -t180 * t193 + t181 * t209, -t167, 0, 0, 0; 0 (t192 * t206 + t207) * t191, t173, 0, 0, 0;];
JR_rot  = t1;
