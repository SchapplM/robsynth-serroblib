% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:23
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRPR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:23:30
% EndTime: 2019-02-22 09:23:30
% DurationCPUTime: 0.14s
% Computational Cost: add. (205->44), mult. (516->91), div. (0->0), fcn. (712->14), ass. (0->47)
t217 = cos(qJ(3));
t191 = pkin(13) + qJ(6);
t189 = sin(t191);
t202 = cos(qJ(4));
t216 = t189 * t202;
t190 = cos(t191);
t215 = t190 * t202;
t193 = sin(pkin(11));
t199 = cos(pkin(6));
t214 = t193 * t199;
t194 = sin(pkin(7));
t213 = t194 * t199;
t195 = sin(pkin(6));
t212 = t195 * t194;
t198 = cos(pkin(7));
t211 = t195 * t198;
t196 = cos(pkin(12));
t210 = t196 * t198;
t197 = cos(pkin(11));
t209 = t197 * t199;
t208 = t195 * t217;
t207 = t194 * t208;
t192 = sin(pkin(12));
t206 = -t192 * t193 + t196 * t209;
t205 = t192 * t197 + t196 * t214;
t204 = t206 * t198;
t203 = t205 * t198;
t201 = sin(qJ(3));
t200 = sin(qJ(4));
t185 = -t192 * t214 + t196 * t197;
t184 = t192 * t209 + t193 * t196;
t183 = -t196 * t212 + t198 * t199;
t180 = t193 * t211 + t194 * t205;
t179 = -t194 * t206 - t197 * t211;
t178 = t201 * t213 + (t192 * t217 + t201 * t210) * t195;
t177 = t192 * t195 * t201 - t208 * t210 - t213 * t217;
t176 = t178 * t202 + t183 * t200;
t175 = -t178 * t200 + t183 * t202;
t174 = t185 * t217 + (t193 * t212 - t203) * t201;
t173 = t185 * t201 - t193 * t207 + t203 * t217;
t172 = t184 * t217 + (-t197 * t212 + t204) * t201;
t171 = t184 * t201 + t197 * t207 - t204 * t217;
t170 = t174 * t202 + t180 * t200;
t169 = -t174 * t200 + t180 * t202;
t168 = t172 * t202 + t179 * t200;
t167 = -t172 * t200 + t179 * t202;
t1 = [0, 0, -t173 * t215 + t174 * t189, t169 * t190, 0, -t170 * t189 + t173 * t190; 0, 0, -t171 * t215 + t172 * t189, t167 * t190, 0, -t168 * t189 + t171 * t190; 0, 0, -t177 * t215 + t178 * t189, t175 * t190, 0, -t176 * t189 + t177 * t190; 0, 0, t173 * t216 + t174 * t190, -t169 * t189, 0, -t170 * t190 - t173 * t189; 0, 0, t171 * t216 + t172 * t190, -t167 * t189, 0, -t168 * t190 - t171 * t189; 0, 0, t177 * t216 + t178 * t190, -t175 * t189, 0, -t176 * t190 - t177 * t189; 0, 0, -t173 * t200, t170, 0, 0; 0, 0, -t171 * t200, t168, 0, 0; 0, 0, -t177 * t200, t176, 0, 0;];
JR_rot  = t1;
