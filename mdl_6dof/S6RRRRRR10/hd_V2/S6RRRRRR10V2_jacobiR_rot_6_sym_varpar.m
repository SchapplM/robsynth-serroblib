% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR10V2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiR_rot_6_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:32
% EndTime: 2019-04-11 14:56:32
% DurationCPUTime: 0.17s
% Computational Cost: add. (220->39), mult. (341->83), div. (0->0), fcn. (490->10), ass. (0->54)
t198 = qJ(2) + qJ(3);
t197 = cos(t198);
t201 = sin(qJ(4));
t206 = cos(qJ(1));
t211 = t206 * t201;
t202 = sin(qJ(1));
t205 = cos(qJ(4));
t214 = t202 * t205;
t189 = t197 * t214 - t211;
t204 = cos(qJ(5));
t196 = sin(t198);
t200 = sin(qJ(5));
t222 = t196 * t200;
t176 = t189 * t204 + t202 * t222;
t210 = t206 * t205;
t215 = t202 * t201;
t188 = t197 * t215 + t210;
t199 = sin(qJ(6));
t203 = cos(qJ(6));
t226 = t176 * t199 - t188 * t203;
t225 = -t176 * t203 - t188 * t199;
t221 = t196 * t204;
t220 = t196 * t206;
t219 = t199 * t201;
t218 = t199 * t204;
t217 = t200 * t205;
t216 = t201 * t203;
t213 = t203 * t204;
t212 = t204 * t205;
t209 = t196 * t219;
t208 = t196 * t216;
t207 = t196 * t211;
t175 = -t189 * t200 + t202 * t221;
t185 = t196 * t212 - t197 * t200;
t184 = -t196 * t217 - t197 * t204;
t191 = t197 * t210 + t215;
t190 = t197 * t211 - t214;
t187 = t197 * t212 + t222;
t186 = t197 * t217 - t221;
t183 = t185 * t206;
t182 = t184 * t206;
t181 = t185 * t202;
t180 = t184 * t202;
t179 = t191 * t204 + t200 * t220;
t178 = t191 * t200 - t204 * t220;
t174 = t187 * t203 + t197 * t219;
t173 = -t187 * t199 + t197 * t216;
t172 = -t183 * t203 - t199 * t207;
t171 = t183 * t199 - t203 * t207;
t170 = -t181 * t203 - t202 * t209;
t169 = t181 * t199 - t202 * t208;
t168 = t179 * t203 + t190 * t199;
t167 = -t179 * t199 + t190 * t203;
t1 = [t225, t172, t172, -t190 * t213 + t191 * t199, -t178 * t203, t167; t168, t170, t170, -t188 * t213 + t189 * t199, t175 * t203, -t226; 0, t174, t174 (t199 * t205 - t201 * t213) * t196, t184 * t203, -t185 * t199 + t208; t226, t171, t171, t190 * t218 + t191 * t203, t178 * t199, -t168; t167, t169, t169, t188 * t218 + t189 * t203, -t175 * t199, t225; 0, t173, t173 (t201 * t218 + t203 * t205) * t196, -t184 * t199, -t185 * t203 - t209; t175, t182, t182, -t190 * t200, t179, 0; t178, t180, t180, -t188 * t200, t176, 0; 0, t186, t186, -t201 * t222, t185, 0;];
JR_rot  = t1;
