% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:48
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR14_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiR_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:47:58
% EndTime: 2019-02-22 11:47:59
% DurationCPUTime: 0.25s
% Computational Cost: add. (180->48), mult. (545->107), div. (0->0), fcn. (746->14), ass. (0->57)
t199 = sin(qJ(2));
t200 = sin(qJ(1));
t202 = cos(qJ(2));
t203 = cos(qJ(1));
t228 = cos(pkin(6));
t214 = t203 * t228;
t186 = t199 * t214 + t200 * t202;
t191 = sin(pkin(14));
t195 = cos(pkin(14));
t185 = t200 * t199 - t202 * t214;
t193 = sin(pkin(7));
t197 = cos(pkin(7));
t194 = sin(pkin(6));
t220 = t194 * t203;
t206 = t185 * t197 + t193 * t220;
t170 = -t186 * t195 + t206 * t191;
t198 = sin(qJ(4));
t201 = cos(qJ(4));
t169 = t186 * t191 + t206 * t195;
t179 = -t185 * t193 + t197 * t220;
t192 = sin(pkin(8));
t196 = cos(pkin(8));
t211 = t169 * t196 + t179 * t192;
t234 = t170 * t201 + t211 * t198;
t233 = -t170 * t198 + t211 * t201;
t225 = t191 * t197;
t224 = t192 * t193;
t223 = t193 * t196;
t222 = t194 * t199;
t221 = t194 * t200;
t219 = t195 * t197;
t218 = t197 * t199;
t217 = t197 * t202;
t216 = t193 * t222;
t215 = t200 * t228;
t213 = t228 * t193;
t188 = -t199 * t215 + t203 * t202;
t187 = -t203 * t199 - t202 * t215;
t205 = t187 * t197 + t193 * t221;
t171 = -t188 * t191 + t205 * t195;
t181 = -t187 * t193 + t197 * t221;
t210 = t171 * t196 + t181 * t192;
t209 = (t195 * t213 + (-t191 * t199 + t195 * t217) * t194) * t196 + (-t194 * t202 * t193 + t228 * t197) * t192;
t173 = t185 * t191 - t186 * t219;
t208 = t173 * t196 + t186 * t224;
t175 = -t187 * t191 - t188 * t219;
t207 = t175 * t196 + t188 * t224;
t182 = (-t191 * t202 - t195 * t218) * t194;
t204 = t182 * t196 + t192 * t216;
t183 = (-t191 * t218 + t195 * t202) * t194;
t178 = t195 * t222 + (t194 * t217 + t213) * t191;
t176 = t187 * t195 - t188 * t225;
t174 = -t185 * t195 - t186 * t225;
t172 = t188 * t195 + t205 * t191;
t166 = t172 * t201 + t210 * t198;
t165 = -t172 * t198 + t210 * t201;
t1 = [t234, t176 * t201 + t207 * t198, 0, t165, 0, 0; t166, t174 * t201 + t208 * t198, 0, -t233, 0, 0; 0, t183 * t201 + t204 * t198, 0, -t178 * t198 + t209 * t201, 0, 0; t233, -t176 * t198 + t207 * t201, 0, -t166, 0, 0; t165, -t174 * t198 + t208 * t201, 0, t234, 0, 0; 0, -t183 * t198 + t204 * t201, 0, -t178 * t201 - t209 * t198, 0, 0; -t169 * t192 + t179 * t196, -t175 * t192 + t188 * t223, 0, 0, 0, 0; -t171 * t192 + t181 * t196, -t173 * t192 + t186 * t223, 0, 0, 0, 0; 0, -t182 * t192 + t196 * t216, 0, 0, 0, 0;];
JR_rot  = t1;
