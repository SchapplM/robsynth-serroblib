% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRP6_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:03
% EndTime: 2019-02-26 20:18:03
% DurationCPUTime: 0.06s
% Computational Cost: add. (34->21), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->29)
t197 = sin(pkin(12));
t199 = sin(pkin(6));
t217 = t197 * t199;
t198 = sin(pkin(7));
t216 = t198 * t199;
t202 = cos(pkin(6));
t215 = t198 * t202;
t200 = cos(pkin(12));
t214 = t200 * t199;
t201 = cos(pkin(7));
t208 = cos(qJ(2));
t213 = t201 * t208;
t205 = sin(qJ(2));
t212 = t202 * t205;
t211 = t202 * t208;
t193 = -t197 * t205 + t200 * t211;
t210 = -t193 * t201 + t198 * t214;
t195 = -t197 * t211 - t200 * t205;
t209 = t195 * t201 + t197 * t216;
t207 = cos(qJ(3));
t206 = cos(qJ(4));
t204 = sin(qJ(3));
t203 = sin(qJ(4));
t196 = -t197 * t212 + t200 * t208;
t194 = t197 * t208 + t200 * t212;
t192 = t202 * t201 - t208 * t216;
t191 = -t195 * t198 + t201 * t217;
t190 = -t193 * t198 - t201 * t214;
t1 = [0, t217, t191, t196 * t204 - t209 * t207 (t196 * t207 + t209 * t204) * t203 - t191 * t206, 0; 0, -t214, t190, t194 * t204 + t210 * t207 (t194 * t207 - t210 * t204) * t203 - t190 * t206, 0; 0, t202, t192, -t207 * t215 + (t204 * t205 - t207 * t213) * t199 (t204 * t215 + (t204 * t213 + t205 * t207) * t199) * t203 - t192 * t206, 0;];
Jg_rot  = t1;
