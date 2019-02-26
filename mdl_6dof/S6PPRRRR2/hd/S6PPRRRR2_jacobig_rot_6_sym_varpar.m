% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRR2_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:16
% EndTime: 2019-02-26 19:43:16
% DurationCPUTime: 0.06s
% Computational Cost: add. (51->20), mult. (150->46), div. (0->0), fcn. (213->12), ass. (0->31)
t157 = sin(pkin(12));
t163 = cos(pkin(6));
t175 = t157 * t163;
t158 = sin(pkin(7));
t159 = sin(pkin(6));
t174 = t158 * t159;
t173 = t158 * t163;
t162 = cos(pkin(7));
t172 = t159 * t162;
t160 = cos(pkin(13));
t171 = t160 * t162;
t161 = cos(pkin(12));
t170 = t161 * t163;
t156 = sin(pkin(13));
t152 = -t157 * t156 + t160 * t170;
t169 = -t152 * t162 + t161 * t174;
t154 = -t161 * t156 - t160 * t175;
t168 = t154 * t162 + t157 * t174;
t167 = cos(qJ(3));
t166 = cos(qJ(4));
t165 = sin(qJ(3));
t164 = sin(qJ(4));
t155 = -t156 * t175 + t161 * t160;
t153 = t156 * t170 + t157 * t160;
t151 = -t160 * t174 + t163 * t162;
t150 = -t154 * t158 + t157 * t172;
t149 = -t152 * t158 - t161 * t172;
t148 = (t165 * t173 + (t156 * t167 + t165 * t171) * t159) * t164 - t151 * t166;
t147 = (t155 * t167 + t168 * t165) * t164 - t150 * t166;
t146 = (t153 * t167 - t169 * t165) * t164 - t149 * t166;
t1 = [0, 0, t150, t155 * t165 - t168 * t167, t147, t147; 0, 0, t149, t153 * t165 + t169 * t167, t146, t146; 0, 0, t151, -t167 * t173 + (t156 * t165 - t167 * t171) * t159, t148, t148;];
Jg_rot  = t1;
