% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRPRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:48
% EndTime: 2019-02-26 19:39:49
% DurationCPUTime: 0.19s
% Computational Cost: add. (283->42), mult. (804->95), div. (0->0), fcn. (1108->16), ass. (0->49)
t203 = sin(pkin(11));
t205 = sin(pkin(6));
t227 = t203 * t205;
t210 = cos(pkin(6));
t226 = t203 * t210;
t208 = cos(pkin(11));
t225 = t205 * t208;
t209 = cos(pkin(7));
t224 = t205 * t209;
t223 = t208 * t210;
t211 = sin(qJ(6));
t215 = cos(qJ(5));
t222 = t211 * t215;
t214 = cos(qJ(6));
t221 = t214 * t215;
t201 = sin(pkin(13));
t206 = cos(pkin(13));
t213 = sin(qJ(3));
t216 = cos(qJ(3));
t220 = t216 * t201 + t213 * t206;
t198 = t213 * t201 - t216 * t206;
t204 = sin(pkin(7));
t190 = t220 * t204;
t192 = t220 * t209;
t202 = sin(pkin(12));
t207 = cos(pkin(12));
t194 = -t203 * t202 + t207 * t223;
t195 = t202 * t223 + t203 * t207;
t219 = -t190 * t225 + t194 * t192 - t195 * t198;
t196 = -t208 * t202 - t207 * t226;
t197 = -t202 * t226 + t208 * t207;
t218 = t190 * t227 + t196 * t192 - t197 * t198;
t217 = t210 * t190 + (t192 * t207 - t198 * t202) * t205;
t212 = sin(qJ(5));
t193 = -t205 * t207 * t204 + t210 * t209;
t191 = t198 * t209;
t189 = t198 * t204;
t187 = -t196 * t204 + t203 * t224;
t186 = -t194 * t204 - t208 * t224;
t184 = -t210 * t189 + (-t191 * t207 - t202 * t220) * t205;
t182 = t193 * t212 + t215 * t217;
t181 = t193 * t215 - t212 * t217;
t179 = -t189 * t227 - t196 * t191 - t197 * t220;
t176 = t189 * t225 - t194 * t191 - t195 * t220;
t174 = t187 * t212 + t215 * t218;
t173 = t187 * t215 - t212 * t218;
t172 = t186 * t212 + t215 * t219;
t171 = t186 * t215 - t212 * t219;
t1 = [0, 0, t179 * t221 + t211 * t218, 0, t173 * t214, -t174 * t211 - t179 * t214; 0, 0, t176 * t221 + t211 * t219, 0, t171 * t214, -t172 * t211 - t176 * t214; 0, 0, t184 * t221 + t211 * t217, 0, t181 * t214, -t182 * t211 - t184 * t214; 0, 0, -t179 * t222 + t214 * t218, 0, -t173 * t211, -t174 * t214 + t179 * t211; 0, 0, -t176 * t222 + t214 * t219, 0, -t171 * t211, -t172 * t214 + t176 * t211; 0, 0, -t184 * t222 + t214 * t217, 0, -t181 * t211, -t182 * t214 + t184 * t211; 0, 0, t179 * t212, 0, t174, 0; 0, 0, t176 * t212, 0, t172, 0; 0, 0, t184 * t212, 0, t182, 0;];
JR_rot  = t1;
