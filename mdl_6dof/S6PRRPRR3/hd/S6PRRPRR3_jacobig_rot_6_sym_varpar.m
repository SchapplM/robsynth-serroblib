% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:14
% EndTime: 2019-02-26 20:05:14
% DurationCPUTime: 0.10s
% Computational Cost: add. (52->25), mult. (148->53), div. (0->0), fcn. (211->14), ass. (0->32)
t181 = sin(pkin(12));
t183 = sin(pkin(6));
t198 = t181 * t183;
t185 = cos(pkin(12));
t197 = t185 * t183;
t187 = cos(pkin(6));
t190 = sin(qJ(2));
t196 = t187 * t190;
t193 = cos(qJ(2));
t195 = t187 * t193;
t180 = sin(pkin(13));
t184 = cos(pkin(13));
t189 = sin(qJ(3));
t192 = cos(qJ(3));
t194 = t192 * t180 + t189 * t184;
t179 = -t189 * t180 + t192 * t184;
t191 = cos(qJ(5));
t188 = sin(qJ(5));
t186 = cos(pkin(7));
t182 = sin(pkin(7));
t177 = -t181 * t196 + t185 * t193;
t176 = -t181 * t195 - t185 * t190;
t175 = t181 * t193 + t185 * t196;
t174 = -t181 * t190 + t185 * t195;
t173 = -t183 * t193 * t182 + t187 * t186;
t172 = t194 * t186;
t171 = t179 * t186;
t170 = t194 * t182;
t169 = t179 * t182;
t168 = -t176 * t182 + t186 * t198;
t167 = -t174 * t182 - t186 * t197;
t1 = [0, t198, t168, 0, -t169 * t198 - t176 * t171 + t177 * t194 (t170 * t198 + t176 * t172 + t177 * t179) * t188 - t168 * t191; 0, -t197, t167, 0, t169 * t197 - t174 * t171 + t175 * t194 (-t170 * t197 + t174 * t172 + t175 * t179) * t188 - t167 * t191; 0, t187, t173, 0, -t187 * t169 + (-t171 * t193 + t190 * t194) * t183 (t187 * t170 + (t172 * t193 + t179 * t190) * t183) * t188 - t173 * t191;];
Jg_rot  = t1;
