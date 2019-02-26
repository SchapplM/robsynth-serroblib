% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR8
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
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR8_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:33
% EndTime: 2019-02-26 22:51:33
% DurationCPUTime: 0.04s
% Computational Cost: add. (26->14), mult. (79->31), div. (0->0), fcn. (116->10), ass. (0->22)
t158 = sin(pkin(6));
t163 = sin(qJ(1));
t172 = t163 * t158;
t162 = sin(qJ(2));
t171 = t163 * t162;
t165 = cos(qJ(2));
t170 = t163 * t165;
t166 = cos(qJ(1));
t169 = t166 * t158;
t168 = t166 * t162;
t167 = t166 * t165;
t164 = cos(qJ(3));
t161 = sin(qJ(3));
t160 = cos(pkin(6));
t159 = cos(pkin(7));
t157 = sin(pkin(7));
t156 = -t160 * t170 - t168;
t155 = t160 * t167 - t171;
t154 = -t160 * t157 * t164 + (-t159 * t164 * t165 + t161 * t162) * t158;
t153 = (-t160 * t171 + t167) * t161 + (-t156 * t159 - t157 * t172) * t164;
t152 = (t160 * t168 + t170) * t161 + (-t155 * t159 + t157 * t169) * t164;
t1 = [0, t172, -t156 * t157 + t159 * t172, t153, t153, 0; 0, -t169, -t155 * t157 - t159 * t169, t152, t152, 0; 1, t160, -t158 * t165 * t157 + t160 * t159, t154, t154, 0;];
Jg_rot  = t1;
