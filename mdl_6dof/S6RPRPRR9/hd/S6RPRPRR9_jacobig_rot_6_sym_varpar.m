% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRPRR9_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:24
% EndTime: 2019-02-26 20:53:24
% DurationCPUTime: 0.09s
% Computational Cost: add. (51->24), mult. (146->51), div. (0->0), fcn. (206->14), ass. (0->34)
t182 = sin(pkin(6));
t189 = sin(qJ(1));
t199 = t182 * t189;
t192 = cos(qJ(1));
t198 = t182 * t192;
t180 = sin(pkin(12));
t197 = t189 * t180;
t184 = cos(pkin(12));
t196 = t189 * t184;
t195 = t192 * t180;
t194 = t192 * t184;
t179 = sin(pkin(13));
t183 = cos(pkin(13));
t188 = sin(qJ(3));
t191 = cos(qJ(3));
t193 = t191 * t179 + t188 * t183;
t178 = -t188 * t179 + t191 * t183;
t190 = cos(qJ(5));
t187 = sin(qJ(5));
t186 = cos(pkin(6));
t185 = cos(pkin(7));
t181 = sin(pkin(7));
t176 = -t186 * t197 + t194;
t175 = -t186 * t196 - t195;
t174 = t186 * t195 + t196;
t173 = t186 * t194 - t197;
t172 = -t182 * t184 * t181 + t186 * t185;
t171 = t193 * t185;
t170 = t178 * t185;
t169 = t193 * t181;
t168 = t178 * t181;
t167 = -t175 * t181 + t185 * t199;
t166 = -t173 * t181 - t185 * t198;
t1 = [0, 0, t167, 0, -t168 * t199 - t175 * t170 + t176 * t193 (t169 * t199 + t175 * t171 + t176 * t178) * t187 - t167 * t190; 0, 0, t166, 0, t168 * t198 - t173 * t170 + t174 * t193 (-t169 * t198 + t173 * t171 + t174 * t178) * t187 - t166 * t190; 1, 0, t172, 0, -t186 * t168 + (-t170 * t184 + t180 * t193) * t182 (t186 * t169 + (t171 * t184 + t178 * t180) * t182) * t187 - t172 * t190;];
Jg_rot  = t1;
