% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR4_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:46
% EndTime: 2019-02-26 21:55:46
% DurationCPUTime: 0.05s
% Computational Cost: add. (31->11), mult. (70->24), div. (0->0), fcn. (106->10), ass. (0->21)
t159 = sin(pkin(6));
t163 = sin(qJ(1));
t169 = t163 * t159;
t165 = cos(qJ(1));
t168 = t165 * t159;
t158 = sin(pkin(12));
t160 = cos(pkin(12));
t162 = sin(qJ(2));
t164 = cos(qJ(2));
t167 = t164 * t158 + t162 * t160;
t166 = t162 * t158 - t164 * t160;
t161 = cos(pkin(6));
t157 = qJ(4) + qJ(5);
t156 = cos(t157);
t155 = sin(t157);
t152 = t167 * t161;
t151 = t166 * t161;
t150 = t166 * t159;
t149 = -t163 * t151 + t165 * t167;
t148 = t165 * t151 + t163 * t167;
t1 = [0, t169, 0, t149, t149 (-t163 * t152 - t165 * t166) * t155 - t156 * t169; 0, -t168, 0, t148, t148 (t165 * t152 - t163 * t166) * t155 + t156 * t168; 1, t161, 0, t150, t150, t167 * t155 * t159 - t161 * t156;];
Jg_rot  = t1;
