% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:10
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR13_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:10:53
% EndTime: 2019-02-22 12:10:53
% DurationCPUTime: 0.20s
% Computational Cost: add. (174->43), mult. (408->90), div. (0->0), fcn. (571->12), ass. (0->51)
t193 = sin(qJ(2));
t194 = sin(qJ(1));
t196 = cos(qJ(2));
t197 = cos(qJ(1));
t216 = cos(pkin(6));
t199 = t197 * t216;
t180 = t193 * t199 + t194 * t196;
t192 = sin(qJ(3));
t195 = cos(qJ(3));
t179 = t194 * t193 - t196 * t199;
t189 = sin(pkin(7));
t191 = cos(pkin(7));
t190 = sin(pkin(6));
t210 = t190 * t197;
t198 = t179 * t191 + t189 * t210;
t166 = -t180 * t195 + t198 * t192;
t173 = -t179 * t189 + t191 * t210;
t188 = pkin(13) + qJ(5);
t186 = sin(t188);
t187 = cos(t188);
t220 = t166 * t186 - t173 * t187;
t219 = t166 * t187 + t173 * t186;
t214 = t186 * t189;
t213 = t187 * t189;
t212 = t189 * t190;
t211 = t190 * t194;
t209 = t191 * t192;
t208 = t191 * t195;
t207 = t192 * t193;
t206 = t192 * t196;
t205 = t193 * t195;
t204 = t195 * t196;
t203 = t193 * t212;
t202 = t189 * t211;
t201 = t189 * t216;
t200 = t194 * t216;
t164 = -t180 * t192 - t198 * t195;
t182 = -t193 * t200 + t197 * t196;
t181 = -t197 * t193 - t196 * t200;
t178 = t216 * t191 - t196 * t212;
t177 = (-t191 * t207 + t204) * t190;
t175 = -t181 * t189 + t191 * t211;
t172 = t192 * t201 + (t191 * t206 + t205) * t190;
t171 = t195 * t201 + (t191 * t204 - t207) * t190;
t170 = t181 * t195 - t182 * t209;
t169 = -t179 * t195 - t180 * t209;
t168 = t182 * t195 + (t181 * t191 + t202) * t192;
t167 = -t181 * t208 + t182 * t192 - t195 * t202;
t163 = t168 * t187 + t175 * t186;
t162 = -t168 * t186 + t175 * t187;
t1 = [t219, t170 * t187 + t182 * t214, -t167 * t187, 0, t162, 0; t163, t169 * t187 + t180 * t214, t164 * t187, 0, t220, 0; 0, t177 * t187 + t186 * t203, t171 * t187, 0, -t172 * t186 + t178 * t187, 0; -t220, -t170 * t186 + t182 * t213, t167 * t186, 0, -t163, 0; t162, -t169 * t186 + t180 * t213, -t164 * t186, 0, t219, 0; 0, -t177 * t186 + t187 * t203, -t171 * t186, 0, -t172 * t187 - t178 * t186, 0; t164, t181 * t192 + t182 * t208, t168, 0, 0, 0; t167, -t179 * t192 + t180 * t208, -t166, 0, 0, 0; 0 (t191 * t205 + t206) * t190, t172, 0, 0, 0;];
JR_rot  = t1;
