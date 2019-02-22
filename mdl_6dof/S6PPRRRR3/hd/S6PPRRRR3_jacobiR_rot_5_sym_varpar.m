% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:27
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiR_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:27:22
% EndTime: 2019-02-22 09:27:22
% DurationCPUTime: 0.14s
% Computational Cost: add. (295->51), mult. (875->112), div. (0->0), fcn. (1191->16), ass. (0->57)
t205 = sin(pkin(13));
t213 = cos(pkin(6));
t234 = t205 * t213;
t206 = sin(pkin(8));
t214 = sin(qJ(5));
t233 = t206 * t214;
t217 = cos(qJ(5));
t232 = t206 * t217;
t207 = sin(pkin(7));
t208 = sin(pkin(6));
t231 = t207 * t208;
t230 = t207 * t213;
t212 = cos(pkin(7));
t229 = t208 * t212;
t210 = cos(pkin(13));
t228 = t210 * t213;
t211 = cos(pkin(8));
t215 = sin(qJ(4));
t227 = t211 * t215;
t218 = cos(qJ(4));
t226 = t211 * t218;
t216 = sin(qJ(3));
t225 = t212 * t216;
t224 = t210 * t231;
t204 = sin(pkin(14));
t209 = cos(pkin(14));
t199 = -t205 * t204 + t209 * t228;
t200 = t204 * t228 + t205 * t209;
t219 = cos(qJ(3));
t189 = -t200 * t216 + (t199 * t212 - t224) * t219;
t196 = -t199 * t207 - t210 * t229;
t223 = t189 * t211 + t196 * t206;
t202 = -t204 * t234 + t210 * t209;
t201 = -t210 * t204 - t209 * t234;
t220 = t201 * t212 + t205 * t231;
t191 = -t202 * t216 + t220 * t219;
t197 = -t201 * t207 + t205 * t229;
t222 = t191 * t211 + t197 * t206;
t194 = t219 * t230 + (t209 * t212 * t219 - t204 * t216) * t208;
t198 = -t209 * t231 + t213 * t212;
t221 = t194 * t211 + t198 * t206;
t195 = t216 * t230 + (t204 * t219 + t209 * t225) * t208;
t193 = -t194 * t206 + t198 * t211;
t192 = t202 * t219 + t220 * t216;
t190 = t199 * t225 + t200 * t219 - t216 * t224;
t188 = t194 * t218 - t195 * t227;
t187 = -t191 * t206 + t197 * t211;
t186 = -t189 * t206 + t196 * t211;
t185 = t195 * t218 + t221 * t215;
t184 = -t195 * t215 + t221 * t218;
t183 = t191 * t218 - t192 * t227;
t182 = t189 * t218 - t190 * t227;
t181 = t192 * t218 + t222 * t215;
t180 = -t192 * t215 + t222 * t218;
t179 = t190 * t218 + t223 * t215;
t178 = -t190 * t215 + t223 * t218;
t1 = [0, 0, t183 * t217 + t192 * t233, t180 * t217, -t181 * t214 + t187 * t217, 0; 0, 0, t182 * t217 + t190 * t233, t178 * t217, -t179 * t214 + t186 * t217, 0; 0, 0, t188 * t217 + t195 * t233, t184 * t217, -t185 * t214 + t193 * t217, 0; 0, 0, -t183 * t214 + t192 * t232, -t180 * t214, -t181 * t217 - t187 * t214, 0; 0, 0, -t182 * t214 + t190 * t232, -t178 * t214, -t179 * t217 - t186 * t214, 0; 0, 0, -t188 * t214 + t195 * t232, -t184 * t214, -t185 * t217 - t193 * t214, 0; 0, 0, t191 * t215 + t192 * t226, t181, 0, 0; 0, 0, t189 * t215 + t190 * t226, t179, 0, 0; 0, 0, t194 * t215 + t195 * t226, t185, 0, 0;];
JR_rot  = t1;
