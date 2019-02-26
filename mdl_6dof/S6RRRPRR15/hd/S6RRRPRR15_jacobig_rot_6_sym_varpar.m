% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR15_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:16
% EndTime: 2019-02-26 22:24:16
% DurationCPUTime: 0.08s
% Computational Cost: add. (34->21), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->30)
t188 = sin(pkin(7));
t191 = cos(pkin(6));
t209 = t188 * t191;
t190 = cos(pkin(7));
t198 = cos(qJ(2));
t208 = t190 * t198;
t189 = sin(pkin(6));
t195 = sin(qJ(1));
t207 = t195 * t189;
t194 = sin(qJ(2));
t206 = t195 * t194;
t205 = t195 * t198;
t199 = cos(qJ(1));
t204 = t199 * t189;
t203 = t199 * t194;
t202 = t199 * t198;
t184 = t191 * t202 - t206;
t201 = -t184 * t190 + t188 * t204;
t186 = -t191 * t205 - t203;
t200 = t186 * t190 + t188 * t207;
t197 = cos(qJ(3));
t196 = cos(qJ(5));
t193 = sin(qJ(3));
t192 = sin(qJ(5));
t187 = -t191 * t206 + t202;
t185 = t191 * t203 + t205;
t183 = -t189 * t198 * t188 + t191 * t190;
t182 = -t186 * t188 + t190 * t207;
t181 = -t184 * t188 - t190 * t204;
t1 = [0, t207, t182, 0, t187 * t197 + t193 * t200, t182 * t192 - (t187 * t193 - t197 * t200) * t196; 0, -t204, t181, 0, t185 * t197 - t201 * t193, t181 * t192 - (t185 * t193 + t197 * t201) * t196; 1, t191, t183, 0, t193 * t209 + (t193 * t208 + t194 * t197) * t189, t183 * t192 - (-t197 * t209 + (t193 * t194 - t197 * t208) * t189) * t196;];
Jg_rot  = t1;
