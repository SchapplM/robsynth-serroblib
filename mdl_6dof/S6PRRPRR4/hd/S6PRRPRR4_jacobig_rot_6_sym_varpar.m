% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR4_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:52
% EndTime: 2019-02-26 20:05:52
% DurationCPUTime: 0.08s
% Computational Cost: add. (19->17), mult. (52->33), div. (0->0), fcn. (81->10), ass. (0->21)
t152 = sin(pkin(6));
t156 = sin(qJ(3));
t166 = t152 * t156;
t159 = cos(qJ(3));
t165 = t152 * t159;
t160 = cos(qJ(2));
t164 = t152 * t160;
t153 = cos(pkin(11));
t163 = t153 * t152;
t154 = cos(pkin(6));
t157 = sin(qJ(2));
t162 = t154 * t157;
t161 = t154 * t160;
t158 = cos(qJ(5));
t155 = sin(qJ(5));
t151 = sin(pkin(11));
t150 = -t151 * t162 + t153 * t160;
t149 = t151 * t161 + t153 * t157;
t148 = t151 * t160 + t153 * t162;
t147 = t151 * t157 - t153 * t161;
t1 = [0, t151 * t152, t149, 0, -t149 (t150 * t159 + t151 * t166) * t155 - (t150 * t156 - t151 * t165) * t158; 0, -t163, t147, 0, -t147 (t148 * t159 - t156 * t163) * t155 - (t148 * t156 + t159 * t163) * t158; 0, t154, -t164, 0, t164 (t154 * t156 + t157 * t165) * t155 - (-t154 * t159 + t157 * t166) * t158;];
Jg_rot  = t1;
