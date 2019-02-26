% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR11_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:40
% EndTime: 2019-02-26 21:06:40
% DurationCPUTime: 0.16s
% Computational Cost: add. (175->38), mult. (514->82), div. (0->0), fcn. (710->14), ass. (0->45)
t195 = cos(pkin(6));
t189 = sin(pkin(12));
t201 = cos(qJ(1));
t209 = t201 * t189;
t193 = cos(pkin(12));
t198 = sin(qJ(1));
t210 = t198 * t193;
t184 = t195 * t209 + t210;
t197 = sin(qJ(3));
t200 = cos(qJ(3));
t208 = t201 * t193;
t211 = t198 * t189;
t183 = -t195 * t208 + t211;
t190 = sin(pkin(7));
t194 = cos(pkin(7));
t191 = sin(pkin(6));
t214 = t191 * t201;
t206 = t183 * t194 + t190 * t214;
t173 = -t184 * t200 + t206 * t197;
t179 = -t183 * t190 + t194 * t214;
t196 = sin(qJ(4));
t199 = cos(qJ(4));
t165 = t173 * t196 - t179 * t199;
t166 = t173 * t199 + t179 * t196;
t205 = t195 * t210 + t209;
t215 = t191 * t198;
t219 = -t190 * t215 + t205 * t194;
t188 = sin(pkin(13));
t217 = t188 * t199;
t216 = t190 * t195;
t192 = cos(pkin(13));
t213 = t192 * t199;
t212 = t193 * t194;
t203 = -t184 * t197 - t206 * t200;
t202 = t205 * t190 + t194 * t215;
t185 = -t195 * t211 + t208;
t182 = -t191 * t193 * t190 + t195 * t194;
t178 = t197 * t216 + (t189 * t200 + t197 * t212) * t191;
t177 = t200 * t216 + (-t189 * t197 + t200 * t212) * t191;
t175 = t185 * t200 - t219 * t197;
t174 = t185 * t197 + t219 * t200;
t169 = -t178 * t196 + t182 * t199;
t168 = t175 * t199 + t202 * t196;
t167 = t175 * t196 - t202 * t199;
t1 = [t166 * t192 + t188 * t203, 0, -t174 * t213 + t175 * t188, -t167 * t192, 0, 0; t168 * t192 + t174 * t188, 0, -t173 * t188 + t203 * t213, t165 * t192, 0, 0; 0, 0, t177 * t213 + t178 * t188, t169 * t192, 0, 0; -t166 * t188 + t192 * t203, 0, t174 * t217 + t175 * t192, t167 * t188, 0, 0; -t168 * t188 + t174 * t192, 0, -t173 * t192 - t203 * t217, -t165 * t188, 0, 0; 0, 0, -t177 * t217 + t178 * t192, -t169 * t188, 0, 0; t165, 0, -t174 * t196, t168, 0, 0; t167, 0, t203 * t196, -t166, 0, 0; 0, 0, t177 * t196, t178 * t199 + t182 * t196, 0, 0;];
JR_rot  = t1;
