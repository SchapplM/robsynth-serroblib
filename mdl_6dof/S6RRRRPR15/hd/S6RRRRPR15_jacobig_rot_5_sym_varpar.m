% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Jg_rot = S6RRRRPR15_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:57
% EndTime: 2019-02-26 22:38:57
% DurationCPUTime: 0.04s
% Computational Cost: add. (16->14), mult. (48->31), div. (0->0), fcn. (72->10), ass. (0->19)
t165 = sin(pkin(6));
t170 = sin(qJ(1));
t179 = t170 * t165;
t169 = sin(qJ(2));
t178 = t170 * t169;
t172 = cos(qJ(2));
t177 = t170 * t172;
t173 = cos(qJ(1));
t176 = t173 * t165;
t175 = t173 * t169;
t174 = t173 * t172;
t171 = cos(qJ(3));
t168 = sin(qJ(3));
t167 = cos(pkin(6));
t166 = cos(pkin(7));
t164 = sin(pkin(7));
t163 = -t167 * t177 - t175;
t162 = t167 * t174 - t178;
t1 = [0, t179, -t163 * t164 + t166 * t179 (-t167 * t178 + t174) * t168 + (-t163 * t166 - t164 * t179) * t171, 0, 0; 0, -t176, -t162 * t164 - t166 * t176 (t167 * t175 + t177) * t168 + (-t162 * t166 + t164 * t176) * t171, 0, 0; 1, t167, -t165 * t172 * t164 + t167 * t166, -t167 * t164 * t171 + (-t166 * t171 * t172 + t168 * t169) * t165, 0, 0;];
Jg_rot  = t1;
