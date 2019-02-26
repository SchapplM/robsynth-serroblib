% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR6_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobig_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:22:02
% EndTime: 2019-02-26 20:22:02
% DurationCPUTime: 0.06s
% Computational Cost: add. (55->27), mult. (161->59), div. (0->0), fcn. (227->14), ass. (0->34)
t201 = sin(pkin(14));
t204 = sin(pkin(6));
t223 = t201 * t204;
t203 = sin(pkin(7));
t222 = t203 * t204;
t208 = cos(pkin(6));
t221 = t203 * t208;
t205 = cos(pkin(14));
t220 = t205 * t204;
t207 = cos(pkin(7));
t214 = cos(qJ(2));
t219 = t207 * t214;
t211 = sin(qJ(2));
t218 = t208 * t211;
t217 = t208 * t214;
t197 = -t201 * t211 + t205 * t217;
t216 = t197 * t207 - t203 * t220;
t199 = -t201 * t217 - t205 * t211;
t215 = t199 * t207 + t201 * t222;
t213 = cos(qJ(3));
t212 = cos(qJ(4));
t210 = sin(qJ(3));
t209 = sin(qJ(4));
t206 = cos(pkin(8));
t202 = sin(pkin(8));
t200 = -t201 * t218 + t205 * t214;
t198 = t201 * t214 + t205 * t218;
t196 = t208 * t207 - t214 * t222;
t195 = -t199 * t203 + t207 * t223;
t194 = -t197 * t203 - t207 * t220;
t193 = t213 * t221 + (-t210 * t211 + t213 * t219) * t204;
t192 = -t200 * t210 + t215 * t213;
t191 = -t198 * t210 + t216 * t213;
t1 = [0, t223, t195, -t192 * t202 + t195 * t206 (t200 * t213 + t215 * t210) * t209 + (-t192 * t206 - t195 * t202) * t212, 0; 0, -t220, t194, -t191 * t202 + t194 * t206 (t198 * t213 + t216 * t210) * t209 + (-t191 * t206 - t194 * t202) * t212, 0; 0, t208, t196, -t193 * t202 + t196 * t206 (t210 * t221 + (t210 * t219 + t211 * t213) * t204) * t209 + (-t193 * t206 - t196 * t202) * t212, 0;];
Jg_rot  = t1;
