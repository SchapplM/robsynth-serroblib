% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPR12_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:37:00
% EndTime: 2019-02-26 22:37:00
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->22), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->31)
t199 = sin(pkin(7));
t202 = cos(pkin(6));
t218 = t199 * t202;
t201 = cos(pkin(7));
t207 = cos(qJ(2));
t217 = t201 * t207;
t200 = sin(pkin(6));
t205 = sin(qJ(1));
t216 = t205 * t200;
t204 = sin(qJ(2));
t215 = t205 * t204;
t214 = t205 * t207;
t208 = cos(qJ(1));
t213 = t208 * t200;
t212 = t208 * t204;
t211 = t208 * t207;
t192 = t202 * t211 - t215;
t210 = -t192 * t201 + t199 * t213;
t194 = -t202 * t214 - t212;
t209 = t194 * t201 + t199 * t216;
t206 = cos(qJ(3));
t203 = sin(qJ(3));
t198 = qJ(4) + pkin(13);
t197 = cos(t198);
t196 = sin(t198);
t195 = -t202 * t215 + t211;
t193 = t202 * t212 + t214;
t191 = -t200 * t207 * t199 + t202 * t201;
t190 = -t194 * t199 + t201 * t216;
t189 = -t192 * t199 - t201 * t213;
t1 = [0, t216, t190, t195 * t203 - t209 * t206, 0 (t195 * t206 + t209 * t203) * t196 - t190 * t197; 0, -t213, t189, t193 * t203 + t210 * t206, 0 (t193 * t206 - t210 * t203) * t196 - t189 * t197; 1, t202, t191, -t206 * t218 + (t203 * t204 - t206 * t217) * t200, 0 (t203 * t218 + (t203 * t217 + t204 * t206) * t200) * t196 - t191 * t197;];
Jg_rot  = t1;
