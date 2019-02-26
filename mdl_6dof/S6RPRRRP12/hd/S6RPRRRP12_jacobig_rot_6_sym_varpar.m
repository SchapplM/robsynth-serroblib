% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRRP12_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:18
% EndTime: 2019-02-26 21:14:18
% DurationCPUTime: 0.09s
% Computational Cost: add. (33->20), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->30)
t193 = sin(pkin(7));
t197 = cos(pkin(6));
t213 = t193 * t197;
t194 = sin(pkin(6));
t200 = sin(qJ(1));
t212 = t194 * t200;
t203 = cos(qJ(1));
t211 = t194 * t203;
t195 = cos(pkin(12));
t196 = cos(pkin(7));
t210 = t195 * t196;
t192 = sin(pkin(12));
t209 = t200 * t192;
t208 = t200 * t195;
t207 = t203 * t192;
t206 = t203 * t195;
t188 = t197 * t206 - t209;
t205 = -t188 * t196 + t193 * t211;
t190 = -t197 * t208 - t207;
t204 = t190 * t196 + t193 * t212;
t202 = cos(qJ(3));
t201 = cos(qJ(4));
t199 = sin(qJ(3));
t198 = sin(qJ(4));
t191 = -t197 * t209 + t206;
t189 = t197 * t207 + t208;
t187 = -t194 * t195 * t193 + t197 * t196;
t186 = -t190 * t193 + t196 * t212;
t185 = -t188 * t193 - t196 * t211;
t1 = [0, 0, t186, t191 * t199 - t204 * t202 (t191 * t202 + t204 * t199) * t198 - t186 * t201, 0; 0, 0, t185, t189 * t199 + t205 * t202 (t189 * t202 - t205 * t199) * t198 - t185 * t201, 0; 1, 0, t187, -t202 * t213 + (t192 * t199 - t202 * t210) * t194 (t199 * t213 + (t192 * t202 + t199 * t210) * t194) * t198 - t187 * t201, 0;];
Jg_rot  = t1;
