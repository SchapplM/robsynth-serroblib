% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR9
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
% Datum: 2019-02-22 12:08
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR9_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:08:15
% EndTime: 2019-02-22 12:08:15
% DurationCPUTime: 0.19s
% Computational Cost: add. (217->46), mult. (621->94), div. (0->0), fcn. (865->14), ass. (0->49)
t198 = sin(pkin(7));
t197 = sin(pkin(13));
t202 = sin(qJ(3));
t218 = cos(pkin(13));
t220 = cos(qJ(3));
t209 = t220 * t197 + t202 * t218;
t184 = t209 * t198;
t200 = cos(pkin(7));
t186 = t209 * t200;
t203 = sin(qJ(2));
t204 = sin(qJ(1));
t206 = cos(qJ(2));
t207 = cos(qJ(1));
t219 = cos(pkin(6));
t210 = t207 * t219;
t188 = t204 * t203 - t206 * t210;
t189 = t203 * t210 + t204 * t206;
t208 = -t202 * t197 + t220 * t218;
t199 = sin(pkin(6));
t213 = t199 * t207;
t172 = t184 * t213 + t188 * t186 - t189 * t208;
t180 = -t188 * t198 + t200 * t213;
t201 = sin(qJ(5));
t205 = cos(qJ(5));
t222 = t172 * t205 + t180 * t201;
t221 = t172 * t201 - t180 * t205;
t217 = t198 * t199;
t216 = t198 * t201;
t215 = t198 * t205;
t214 = t199 * t204;
t212 = t203 * t217;
t211 = t204 * t219;
t183 = t208 * t198;
t185 = t208 * t200;
t170 = -t183 * t213 - t188 * t185 - t189 * t209;
t190 = -t207 * t203 - t206 * t211;
t191 = -t203 * t211 + t207 * t206;
t174 = t184 * t214 + t190 * t186 + t191 * t208;
t176 = t219 * t184 + (t186 * t206 + t203 * t208) * t199;
t187 = t219 * t200 - t206 * t217;
t182 = -t190 * t198 + t200 * t214;
t179 = (-t186 * t203 + t206 * t208) * t199;
t178 = -t191 * t186 + t190 * t208;
t177 = -t189 * t186 - t188 * t208;
t175 = t219 * t183 + (t185 * t206 - t203 * t209) * t199;
t173 = t183 * t214 + t190 * t185 - t191 * t209;
t169 = t174 * t205 + t182 * t201;
t168 = -t174 * t201 + t182 * t205;
t1 = [t222, t178 * t205 + t191 * t216, t173 * t205, 0, t168, 0; t169, t177 * t205 + t189 * t216, t170 * t205, 0, t221, 0; 0, t179 * t205 + t201 * t212, t175 * t205, 0, -t176 * t201 + t187 * t205, 0; -t221, -t178 * t201 + t191 * t215, -t173 * t201, 0, -t169, 0; t168, -t177 * t201 + t189 * t215, -t170 * t201, 0, t222, 0; 0, -t179 * t201 + t205 * t212, -t175 * t201, 0, -t176 * t205 - t187 * t201, 0; t170, t191 * t185 + t190 * t209, t174, 0, 0, 0; -t173, t189 * t185 - t188 * t209, -t172, 0, 0, 0; 0 (t185 * t203 + t206 * t209) * t199, t176, 0, 0, 0;];
JR_rot  = t1;
