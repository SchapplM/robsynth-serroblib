% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobiR_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:43
% EndTime: 2019-02-26 19:42:43
% DurationCPUTime: 0.16s
% Computational Cost: add. (289->47), mult. (678->91), div. (0->0), fcn. (937->14), ass. (0->53)
t239 = cos(qJ(3));
t211 = sin(pkin(13));
t212 = sin(pkin(12));
t215 = cos(pkin(13));
t216 = cos(pkin(12));
t218 = cos(pkin(6));
t228 = t216 * t218;
t203 = t211 * t228 + t212 * t215;
t220 = sin(qJ(3));
t217 = cos(pkin(7));
t225 = -t212 * t211 + t215 * t228;
t223 = t225 * t217;
t213 = sin(pkin(7));
t214 = sin(pkin(6));
t231 = t214 * t213;
t193 = t203 * t239 + (-t216 * t231 + t223) * t220;
t230 = t214 * t217;
t198 = -t225 * t213 - t216 * t230;
t210 = qJ(4) + qJ(5);
t208 = sin(t210);
t209 = cos(t210);
t185 = -t193 * t208 + t198 * t209;
t219 = sin(qJ(6));
t238 = t185 * t219;
t233 = t212 * t218;
t204 = -t211 * t233 + t216 * t215;
t224 = t216 * t211 + t215 * t233;
t222 = t224 * t217;
t195 = t204 * t239 + (t212 * t231 - t222) * t220;
t199 = t212 * t230 + t224 * t213;
t187 = -t195 * t208 + t199 * t209;
t237 = t187 * t219;
t229 = t215 * t217;
t232 = t213 * t218;
t197 = t220 * t232 + (t239 * t211 + t220 * t229) * t214;
t202 = -t215 * t231 + t218 * t217;
t190 = -t197 * t208 + t202 * t209;
t236 = t190 * t219;
t235 = t209 * t219;
t221 = cos(qJ(6));
t234 = t209 * t221;
t227 = t214 * t239;
t226 = t213 * t227;
t196 = t214 * t211 * t220 - t227 * t229 - t239 * t232;
t194 = t204 * t220 - t212 * t226 + t239 * t222;
t192 = t203 * t220 + t216 * t226 - t239 * t223;
t191 = t197 * t209 + t202 * t208;
t189 = t190 * t221;
t188 = t195 * t209 + t199 * t208;
t186 = t193 * t209 + t198 * t208;
t184 = t187 * t221;
t183 = t185 * t221;
t1 = [0, 0, -t194 * t234 + t195 * t219, t184, t184, -t188 * t219 + t194 * t221; 0, 0, -t192 * t234 + t193 * t219, t183, t183, -t186 * t219 + t192 * t221; 0, 0, -t196 * t234 + t197 * t219, t189, t189, -t191 * t219 + t196 * t221; 0, 0, t194 * t235 + t195 * t221, -t237, -t237, -t188 * t221 - t194 * t219; 0, 0, t192 * t235 + t193 * t221, -t238, -t238, -t186 * t221 - t192 * t219; 0, 0, t196 * t235 + t197 * t221, -t236, -t236, -t191 * t221 - t196 * t219; 0, 0, -t194 * t208, t188, t188, 0; 0, 0, -t192 * t208, t186, t186, 0; 0, 0, -t196 * t208, t191, t191, 0;];
JR_rot  = t1;
