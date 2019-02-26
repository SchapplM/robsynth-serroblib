% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRRR10_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:50
% EndTime: 2019-02-26 21:19:50
% DurationCPUTime: 0.06s
% Computational Cost: add. (49->21), mult. (129->45), div. (0->0), fcn. (184->12), ass. (0->34)
t192 = sin(pkin(7));
t196 = cos(pkin(6));
t210 = t192 * t196;
t193 = sin(pkin(6));
t198 = sin(qJ(1));
t209 = t193 * t198;
t200 = cos(qJ(1));
t208 = t193 * t200;
t194 = cos(pkin(13));
t195 = cos(pkin(7));
t207 = t194 * t195;
t191 = sin(pkin(13));
t206 = t198 * t191;
t205 = t198 * t194;
t204 = t200 * t191;
t203 = t200 * t194;
t184 = t196 * t203 - t206;
t202 = -t184 * t195 + t192 * t208;
t186 = -t196 * t205 - t204;
t201 = t186 * t195 + t192 * t209;
t199 = cos(qJ(3));
t197 = sin(qJ(3));
t190 = qJ(4) + qJ(5);
t189 = cos(t190);
t188 = sin(t190);
t187 = -t196 * t206 + t203;
t185 = t196 * t204 + t205;
t183 = -t193 * t194 * t192 + t196 * t195;
t182 = -t186 * t192 + t195 * t209;
t181 = -t184 * t192 - t195 * t208;
t180 = -t199 * t210 + (t191 * t197 - t199 * t207) * t193;
t179 = t187 * t197 - t201 * t199;
t178 = t185 * t197 + t202 * t199;
t1 = [0, 0, t182, t179, t179 (t187 * t199 + t201 * t197) * t188 - t182 * t189; 0, 0, t181, t178, t178 (t185 * t199 - t202 * t197) * t188 - t181 * t189; 1, 0, t183, t180, t180 (t197 * t210 + (t191 * t199 + t197 * t207) * t193) * t188 - t183 * t189;];
Jg_rot  = t1;
