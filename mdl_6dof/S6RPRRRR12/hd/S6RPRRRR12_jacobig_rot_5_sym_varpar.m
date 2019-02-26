% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRRR12_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobig_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:10
% EndTime: 2019-02-26 21:21:10
% DurationCPUTime: 0.09s
% Computational Cost: add. (54->26), mult. (159->57), div. (0->0), fcn. (222->14), ass. (0->35)
t198 = sin(pkin(7));
t203 = cos(pkin(6));
t219 = t198 * t203;
t199 = sin(pkin(6));
t206 = sin(qJ(1));
t218 = t199 * t206;
t209 = cos(qJ(1));
t217 = t199 * t209;
t200 = cos(pkin(14));
t202 = cos(pkin(7));
t216 = t200 * t202;
t196 = sin(pkin(14));
t215 = t206 * t196;
t214 = t206 * t200;
t213 = t209 * t196;
t212 = t209 * t200;
t192 = t203 * t212 - t215;
t211 = t192 * t202 - t198 * t217;
t194 = -t203 * t214 - t213;
t210 = t194 * t202 + t198 * t218;
t208 = cos(qJ(3));
t207 = cos(qJ(4));
t205 = sin(qJ(3));
t204 = sin(qJ(4));
t201 = cos(pkin(8));
t197 = sin(pkin(8));
t195 = -t203 * t215 + t212;
t193 = t203 * t213 + t214;
t191 = -t198 * t199 * t200 + t202 * t203;
t190 = -t194 * t198 + t202 * t218;
t189 = -t192 * t198 - t202 * t217;
t188 = t208 * t219 + (-t196 * t205 + t208 * t216) * t199;
t187 = -t195 * t205 + t208 * t210;
t186 = -t193 * t205 + t208 * t211;
t1 = [0, 0, t190, -t187 * t197 + t190 * t201 (t195 * t208 + t205 * t210) * t204 + (-t187 * t201 - t190 * t197) * t207, 0; 0, 0, t189, -t186 * t197 + t189 * t201 (t193 * t208 + t205 * t211) * t204 + (-t186 * t201 - t189 * t197) * t207, 0; 1, 0, t191, -t188 * t197 + t191 * t201 (t205 * t219 + (t196 * t208 + t205 * t216) * t199) * t204 + (-t188 * t201 - t191 * t197) * t207, 0;];
Jg_rot  = t1;
