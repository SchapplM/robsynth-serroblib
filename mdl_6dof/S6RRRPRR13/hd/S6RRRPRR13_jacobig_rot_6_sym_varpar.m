% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR13_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:05
% EndTime: 2019-02-26 22:23:05
% DurationCPUTime: 0.10s
% Computational Cost: add. (40->22), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->31)
t198 = sin(pkin(7));
t201 = cos(pkin(6));
t217 = t198 * t201;
t200 = cos(pkin(7));
t206 = cos(qJ(2));
t216 = t200 * t206;
t203 = sin(qJ(2));
t207 = cos(qJ(1));
t215 = t203 * t207;
t199 = sin(pkin(6));
t204 = sin(qJ(1));
t214 = t204 * t199;
t213 = t204 * t203;
t212 = t204 * t206;
t211 = t206 * t207;
t210 = t207 * t199;
t191 = t201 * t211 - t213;
t209 = -t191 * t200 + t198 * t210;
t193 = -t201 * t212 - t215;
t208 = t193 * t200 + t198 * t214;
t205 = cos(qJ(3));
t202 = sin(qJ(3));
t197 = pkin(13) + qJ(5);
t196 = cos(t197);
t195 = sin(t197);
t194 = -t201 * t213 + t211;
t192 = t201 * t215 + t212;
t190 = -t198 * t199 * t206 + t200 * t201;
t189 = -t193 * t198 + t200 * t214;
t188 = -t191 * t198 - t200 * t210;
t1 = [0, t214, t189, 0, t194 * t202 - t208 * t205 (t194 * t205 + t208 * t202) * t195 - t189 * t196; 0, -t210, t188, 0, t192 * t202 + t209 * t205 (t192 * t205 - t209 * t202) * t195 - t188 * t196; 1, t201, t190, 0, -t205 * t217 + (t202 * t203 - t205 * t216) * t199 (t202 * t217 + (t202 * t216 + t203 * t205) * t199) * t195 - t190 * t196;];
Jg_rot  = t1;
