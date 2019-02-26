% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:41:34
% EndTime: 2019-02-26 19:41:34
% DurationCPUTime: 0.15s
% Computational Cost: add. (175->43), mult. (516->91), div. (0->0), fcn. (712->14), ass. (0->46)
t217 = cos(qJ(3));
t191 = sin(pkin(11));
t197 = cos(pkin(6));
t216 = t191 * t197;
t192 = sin(pkin(7));
t215 = t192 * t197;
t193 = sin(pkin(6));
t214 = t193 * t192;
t196 = cos(pkin(7));
t213 = t193 * t196;
t194 = cos(pkin(12));
t212 = t194 * t196;
t195 = cos(pkin(11));
t211 = t195 * t197;
t198 = sin(qJ(5));
t202 = cos(qJ(4));
t210 = t198 * t202;
t201 = cos(qJ(5));
t209 = t201 * t202;
t208 = t193 * t217;
t207 = t192 * t208;
t190 = sin(pkin(12));
t206 = -t191 * t190 + t194 * t211;
t205 = t195 * t190 + t194 * t216;
t204 = t206 * t196;
t203 = t205 * t196;
t200 = sin(qJ(3));
t199 = sin(qJ(4));
t186 = -t190 * t216 + t195 * t194;
t185 = t190 * t211 + t191 * t194;
t184 = -t194 * t214 + t197 * t196;
t181 = t191 * t213 + t192 * t205;
t180 = -t192 * t206 - t195 * t213;
t179 = t200 * t215 + (t217 * t190 + t200 * t212) * t193;
t178 = t193 * t190 * t200 - t208 * t212 - t217 * t215;
t177 = t179 * t202 + t184 * t199;
t176 = -t179 * t199 + t184 * t202;
t175 = t186 * t217 + (t191 * t214 - t203) * t200;
t174 = t186 * t200 - t191 * t207 + t217 * t203;
t173 = t185 * t217 + (-t195 * t214 + t204) * t200;
t172 = t185 * t200 + t195 * t207 - t217 * t204;
t171 = t175 * t202 + t181 * t199;
t170 = -t175 * t199 + t181 * t202;
t169 = t173 * t202 + t180 * t199;
t168 = -t173 * t199 + t180 * t202;
t1 = [0, 0, -t174 * t209 + t175 * t198, t170 * t201, -t171 * t198 + t174 * t201, 0; 0, 0, -t172 * t209 + t173 * t198, t168 * t201, -t169 * t198 + t172 * t201, 0; 0, 0, -t178 * t209 + t179 * t198, t176 * t201, -t177 * t198 + t178 * t201, 0; 0, 0, t174 * t210 + t175 * t201, -t170 * t198, -t171 * t201 - t174 * t198, 0; 0, 0, t172 * t210 + t173 * t201, -t168 * t198, -t169 * t201 - t172 * t198, 0; 0, 0, t178 * t210 + t179 * t201, -t176 * t198, -t177 * t201 - t178 * t198, 0; 0, 0, -t174 * t199, t171, 0, 0; 0, 0, -t172 * t199, t169, 0, 0; 0, 0, -t178 * t199, t177, 0, 0;];
JR_rot  = t1;
