% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPR15_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:48
% EndTime: 2019-02-26 22:38:48
% DurationCPUTime: 0.05s
% Computational Cost: add. (34->21), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->30)
t183 = sin(pkin(7));
t186 = cos(pkin(6));
t204 = t183 * t186;
t185 = cos(pkin(7));
t193 = cos(qJ(2));
t203 = t185 * t193;
t184 = sin(pkin(6));
t190 = sin(qJ(1));
t202 = t190 * t184;
t189 = sin(qJ(2));
t201 = t190 * t189;
t200 = t190 * t193;
t194 = cos(qJ(1));
t199 = t194 * t184;
t198 = t194 * t189;
t197 = t194 * t193;
t179 = t186 * t197 - t201;
t196 = -t179 * t185 + t183 * t199;
t181 = -t186 * t200 - t198;
t195 = t181 * t185 + t183 * t202;
t192 = cos(qJ(3));
t191 = cos(qJ(4));
t188 = sin(qJ(3));
t187 = sin(qJ(4));
t182 = -t186 * t201 + t197;
t180 = t186 * t198 + t200;
t178 = -t184 * t193 * t183 + t186 * t185;
t177 = -t181 * t183 + t185 * t202;
t176 = -t179 * t183 - t185 * t199;
t1 = [0, t202, t177, t182 * t188 - t195 * t192, 0 (t182 * t192 + t195 * t188) * t191 + t177 * t187; 0, -t199, t176, t180 * t188 + t196 * t192, 0 (t180 * t192 - t196 * t188) * t191 + t176 * t187; 1, t186, t178, -t192 * t204 + (t188 * t189 - t192 * t203) * t184, 0 (t188 * t204 + (t188 * t203 + t189 * t192) * t184) * t191 + t178 * t187;];
Jg_rot  = t1;
