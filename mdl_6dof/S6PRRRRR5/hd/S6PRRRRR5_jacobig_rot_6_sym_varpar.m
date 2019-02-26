% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR5
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
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRR5_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:11
% EndTime: 2019-02-26 20:21:11
% DurationCPUTime: 0.06s
% Computational Cost: add. (52->21), mult. (152->47), div. (0->0), fcn. (218->12), ass. (0->32)
t180 = sin(pkin(13));
t182 = sin(pkin(6));
t200 = t180 * t182;
t181 = sin(pkin(7));
t199 = t181 * t182;
t185 = cos(pkin(6));
t198 = t181 * t185;
t183 = cos(pkin(13));
t197 = t183 * t182;
t184 = cos(pkin(7));
t191 = cos(qJ(2));
t196 = t184 * t191;
t188 = sin(qJ(2));
t195 = t185 * t188;
t194 = t185 * t191;
t176 = -t180 * t188 + t183 * t194;
t193 = -t176 * t184 + t181 * t197;
t178 = -t180 * t194 - t183 * t188;
t192 = t178 * t184 + t180 * t199;
t190 = cos(qJ(3));
t189 = cos(qJ(4));
t187 = sin(qJ(3));
t186 = sin(qJ(4));
t179 = -t180 * t195 + t183 * t191;
t177 = t180 * t191 + t183 * t195;
t175 = t185 * t184 - t191 * t199;
t174 = -t178 * t181 + t184 * t200;
t173 = -t176 * t181 - t184 * t197;
t172 = (t187 * t198 + (t187 * t196 + t188 * t190) * t182) * t186 - t175 * t189;
t171 = (t179 * t190 + t192 * t187) * t186 - t174 * t189;
t170 = (t177 * t190 - t193 * t187) * t186 - t173 * t189;
t1 = [0, t200, t174, t179 * t187 - t192 * t190, t171, t171; 0, -t197, t173, t177 * t187 + t193 * t190, t170, t170; 0, t185, t175, -t190 * t198 + (t187 * t188 - t190 * t196) * t182, t172, t172;];
Jg_rot  = t1;
