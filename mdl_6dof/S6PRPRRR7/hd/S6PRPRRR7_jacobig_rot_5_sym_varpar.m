% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR7_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobig_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:30
% EndTime: 2019-02-26 19:57:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (50->27), mult. (146->59), div. (0->0), fcn. (204->14), ass. (0->34)
t179 = sin(pkin(13));
t182 = sin(pkin(6));
t200 = t179 * t182;
t181 = sin(pkin(7));
t199 = t181 * t182;
t187 = cos(pkin(6));
t198 = t181 * t187;
t184 = cos(pkin(13));
t197 = t184 * t182;
t186 = cos(pkin(7));
t191 = cos(qJ(2));
t196 = t186 * t191;
t189 = sin(qJ(2));
t195 = t187 * t189;
t194 = t187 * t191;
t174 = -t179 * t189 + t184 * t194;
t193 = t174 * t186 - t181 * t197;
t176 = -t179 * t194 - t184 * t189;
t192 = t176 * t186 + t179 * t199;
t190 = cos(qJ(4));
t188 = sin(qJ(4));
t185 = cos(pkin(8));
t183 = cos(pkin(14));
t180 = sin(pkin(8));
t178 = sin(pkin(14));
t177 = -t179 * t195 + t184 * t191;
t175 = t179 * t191 + t184 * t195;
t173 = t187 * t186 - t191 * t199;
t172 = -t176 * t181 + t186 * t200;
t171 = -t174 * t181 - t186 * t197;
t170 = t183 * t198 + (-t178 * t189 + t183 * t196) * t182;
t169 = -t177 * t178 + t192 * t183;
t168 = -t175 * t178 + t193 * t183;
t1 = [0, t200, 0, -t169 * t180 + t172 * t185 (t177 * t183 + t192 * t178) * t188 + (-t169 * t185 - t172 * t180) * t190, 0; 0, -t197, 0, -t168 * t180 + t171 * t185 (t175 * t183 + t193 * t178) * t188 + (-t168 * t185 - t171 * t180) * t190, 0; 0, t187, 0, -t170 * t180 + t173 * t185 (t182 * t189 * t183 + (t182 * t196 + t198) * t178) * t188 + (-t170 * t185 - t173 * t180) * t190, 0;];
Jg_rot  = t1;
