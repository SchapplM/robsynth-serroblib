% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR9_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:07
% EndTime: 2019-02-26 22:52:07
% DurationCPUTime: 0.06s
% Computational Cost: add. (34->21), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->30)
t187 = sin(pkin(7));
t190 = cos(pkin(6));
t208 = t187 * t190;
t189 = cos(pkin(7));
t197 = cos(qJ(2));
t207 = t189 * t197;
t188 = sin(pkin(6));
t194 = sin(qJ(1));
t206 = t194 * t188;
t193 = sin(qJ(2));
t205 = t194 * t193;
t204 = t194 * t197;
t198 = cos(qJ(1));
t203 = t198 * t188;
t202 = t198 * t193;
t201 = t198 * t197;
t183 = t190 * t201 - t205;
t200 = -t183 * t189 + t187 * t203;
t185 = -t190 * t204 - t202;
t199 = t185 * t189 + t187 * t206;
t196 = cos(qJ(3));
t195 = cos(qJ(4));
t192 = sin(qJ(3));
t191 = sin(qJ(4));
t186 = -t190 * t205 + t201;
t184 = t190 * t202 + t204;
t182 = -t188 * t197 * t187 + t190 * t189;
t181 = -t185 * t187 + t189 * t206;
t180 = -t183 * t187 - t189 * t203;
t1 = [0, t206, t181, t186 * t192 - t199 * t196 (t186 * t196 + t199 * t192) * t191 - t181 * t195, 0; 0, -t203, t180, t184 * t192 + t200 * t196 (t184 * t196 - t200 * t192) * t191 - t180 * t195, 0; 1, t190, t182, -t196 * t208 + (t192 * t193 - t196 * t207) * t188 (t192 * t208 + (t192 * t207 + t193 * t196) * t188) * t191 - t182 * t195, 0;];
Jg_rot  = t1;
