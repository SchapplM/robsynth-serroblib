% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:05
% EndTime: 2019-02-26 20:14:05
% DurationCPUTime: 0.13s
% Computational Cost: add. (169->50), mult. (513->117), div. (0->0), fcn. (706->14), ass. (0->58)
t195 = sin(pkin(13));
t206 = cos(qJ(4));
t226 = t195 * t206;
t197 = sin(pkin(7));
t198 = sin(pkin(6));
t225 = t197 * t198;
t202 = cos(pkin(6));
t224 = t197 * t202;
t203 = sin(qJ(4));
t223 = t197 * t203;
t222 = t197 * t206;
t201 = cos(pkin(7));
t221 = t198 * t201;
t199 = cos(pkin(13));
t220 = t199 * t206;
t204 = sin(qJ(3));
t219 = t201 * t204;
t207 = cos(qJ(3));
t218 = t201 * t207;
t205 = sin(qJ(2));
t217 = t202 * t205;
t208 = cos(qJ(2));
t216 = t202 * t208;
t215 = t204 * t205;
t214 = t204 * t208;
t213 = t205 * t207;
t212 = t207 * t208;
t211 = t205 * t225;
t196 = sin(pkin(12));
t200 = cos(pkin(12));
t190 = -t196 * t205 + t200 * t216;
t210 = t190 * t201 - t200 * t225;
t192 = -t196 * t216 - t200 * t205;
t209 = t192 * t201 + t196 * t225;
t193 = -t196 * t217 + t200 * t208;
t191 = t196 * t208 + t200 * t217;
t189 = t202 * t201 - t208 * t225;
t188 = (-t201 * t215 + t212) * t198;
t187 = (t201 * t213 + t214) * t198;
t186 = -t192 * t197 + t196 * t221;
t185 = -t190 * t197 - t200 * t221;
t184 = t204 * t224 + (t201 * t214 + t213) * t198;
t183 = t207 * t224 + (t201 * t212 - t215) * t198;
t182 = t188 * t206 + t203 * t211;
t181 = t192 * t207 - t193 * t219;
t180 = t192 * t204 + t193 * t218;
t179 = t190 * t207 - t191 * t219;
t178 = t190 * t204 + t191 * t218;
t177 = -t184 * t203 + t189 * t206;
t176 = t193 * t207 + t209 * t204;
t175 = -t193 * t204 + t209 * t207;
t174 = t191 * t207 + t210 * t204;
t173 = -t191 * t204 + t210 * t207;
t172 = t181 * t206 + t193 * t223;
t171 = t179 * t206 + t191 * t223;
t170 = -t176 * t203 + t186 * t206;
t169 = -t174 * t203 + t185 * t206;
t1 = [0, t172 * t199 + t180 * t195, t175 * t220 + t176 * t195, t170 * t199, 0, 0; 0, t171 * t199 + t178 * t195, t173 * t220 + t174 * t195, t169 * t199, 0, 0; 0, t182 * t199 + t187 * t195, t183 * t220 + t184 * t195, t177 * t199, 0, 0; 0, -t172 * t195 + t180 * t199, -t175 * t226 + t176 * t199, -t170 * t195, 0, 0; 0, -t171 * t195 + t178 * t199, -t173 * t226 + t174 * t199, -t169 * t195, 0, 0; 0, -t182 * t195 + t187 * t199, -t183 * t226 + t184 * t199, -t177 * t195, 0, 0; 0, t181 * t203 - t193 * t222, t175 * t203, t176 * t206 + t186 * t203, 0, 0; 0, t179 * t203 - t191 * t222, t173 * t203, t174 * t206 + t185 * t203, 0, 0; 0, t188 * t203 - t206 * t211, t183 * t203, t184 * t206 + t189 * t203, 0, 0;];
JR_rot  = t1;
