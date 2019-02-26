% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRP5_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:29
% EndTime: 2019-02-26 20:17:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (34->21), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->29)
t162 = sin(pkin(12));
t164 = sin(pkin(6));
t182 = t162 * t164;
t163 = sin(pkin(7));
t181 = t163 * t164;
t167 = cos(pkin(6));
t180 = t163 * t167;
t165 = cos(pkin(12));
t179 = t165 * t164;
t166 = cos(pkin(7));
t173 = cos(qJ(2));
t178 = t166 * t173;
t170 = sin(qJ(2));
t177 = t167 * t170;
t176 = t167 * t173;
t158 = -t162 * t170 + t165 * t176;
t175 = -t158 * t166 + t163 * t179;
t160 = -t162 * t176 - t165 * t170;
t174 = t160 * t166 + t162 * t181;
t172 = cos(qJ(3));
t171 = cos(qJ(4));
t169 = sin(qJ(3));
t168 = sin(qJ(4));
t161 = -t162 * t177 + t165 * t173;
t159 = t162 * t173 + t165 * t177;
t157 = t167 * t166 - t173 * t181;
t156 = -t160 * t163 + t166 * t182;
t155 = -t158 * t163 - t166 * t179;
t1 = [0, t182, t156, t161 * t169 - t174 * t172 (t161 * t172 + t174 * t169) * t168 - t156 * t171, 0; 0, -t179, t155, t159 * t169 + t175 * t172 (t159 * t172 - t175 * t169) * t168 - t155 * t171, 0; 0, t167, t157, -t172 * t180 + (t169 * t170 - t172 * t178) * t164 (t169 * t180 + (t169 * t178 + t170 * t172) * t164) * t168 - t157 * t171, 0;];
Jg_rot  = t1;
