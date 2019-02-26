% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRR3_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobig_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:56
% EndTime: 2019-02-26 19:43:56
% DurationCPUTime: 0.09s
% Computational Cost: add. (54->26), mult. (159->58), div. (0->0), fcn. (222->14), ass. (0->33)
t157 = sin(pkin(13));
t165 = cos(pkin(6));
t177 = t157 * t165;
t159 = sin(pkin(7));
t160 = sin(pkin(6));
t176 = t159 * t160;
t175 = t159 * t165;
t164 = cos(pkin(7));
t174 = t160 * t164;
t161 = cos(pkin(14));
t173 = t161 * t164;
t162 = cos(pkin(13));
t172 = t162 * t165;
t156 = sin(pkin(14));
t152 = -t157 * t156 + t161 * t172;
t171 = t152 * t164 - t162 * t176;
t154 = -t162 * t156 - t161 * t177;
t170 = t154 * t164 + t157 * t176;
t169 = cos(qJ(3));
t168 = cos(qJ(4));
t167 = sin(qJ(3));
t166 = sin(qJ(4));
t163 = cos(pkin(8));
t158 = sin(pkin(8));
t155 = -t156 * t177 + t162 * t161;
t153 = t156 * t172 + t157 * t161;
t151 = -t161 * t176 + t165 * t164;
t150 = -t154 * t159 + t157 * t174;
t149 = -t152 * t159 - t162 * t174;
t148 = t169 * t175 + (-t156 * t167 + t169 * t173) * t160;
t147 = -t155 * t167 + t170 * t169;
t146 = -t153 * t167 + t171 * t169;
t1 = [0, 0, t150, -t147 * t158 + t150 * t163 (t155 * t169 + t170 * t167) * t166 + (-t147 * t163 - t150 * t158) * t168, 0; 0, 0, t149, -t146 * t158 + t149 * t163 (t153 * t169 + t171 * t167) * t166 + (-t146 * t163 - t149 * t158) * t168, 0; 0, 0, t151, -t148 * t158 + t151 * t163 (t167 * t175 + (t156 * t169 + t167 * t173) * t160) * t166 + (-t148 * t163 - t151 * t158) * t168, 0;];
Jg_rot  = t1;
