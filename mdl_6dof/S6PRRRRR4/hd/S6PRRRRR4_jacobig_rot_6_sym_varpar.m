% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR4_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:20:36
% EndTime: 2019-02-26 20:20:37
% DurationCPUTime: 0.06s
% Computational Cost: add. (50->22), mult. (131->47), div. (0->0), fcn. (189->12), ass. (0->33)
t190 = sin(pkin(13));
t192 = sin(pkin(6));
t208 = t190 * t192;
t191 = sin(pkin(7));
t207 = t191 * t192;
t195 = cos(pkin(6));
t206 = t191 * t195;
t193 = cos(pkin(13));
t205 = t193 * t192;
t194 = cos(pkin(7));
t199 = cos(qJ(2));
t204 = t194 * t199;
t197 = sin(qJ(2));
t203 = t195 * t197;
t202 = t195 * t199;
t183 = -t190 * t197 + t193 * t202;
t201 = -t183 * t194 + t191 * t205;
t185 = -t190 * t202 - t193 * t197;
t200 = t185 * t194 + t190 * t207;
t198 = cos(qJ(3));
t196 = sin(qJ(3));
t189 = qJ(4) + qJ(5);
t188 = cos(t189);
t187 = sin(t189);
t186 = -t190 * t203 + t193 * t199;
t184 = t190 * t199 + t193 * t203;
t182 = t195 * t194 - t199 * t207;
t181 = -t185 * t191 + t194 * t208;
t180 = -t183 * t191 - t194 * t205;
t179 = -t198 * t206 + (t196 * t197 - t198 * t204) * t192;
t178 = t186 * t196 - t200 * t198;
t177 = t184 * t196 + t201 * t198;
t1 = [0, t208, t181, t178, t178 (t186 * t198 + t200 * t196) * t187 - t181 * t188; 0, -t205, t180, t177, t177 (t184 * t198 - t201 * t196) * t187 - t180 * t188; 0, t195, t182, t179, t179 (t196 * t206 + (t196 * t204 + t197 * t198) * t192) * t187 - t182 * t188;];
Jg_rot  = t1;
