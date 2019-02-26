% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRPR1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:19
% EndTime: 2019-02-26 19:40:19
% DurationCPUTime: 0.08s
% Computational Cost: add. (33->20), mult. (98->46), div. (0->0), fcn. (140->12), ass. (0->28)
t148 = sin(pkin(11));
t154 = cos(pkin(6));
t166 = t148 * t154;
t149 = sin(pkin(7));
t150 = sin(pkin(6));
t165 = t149 * t150;
t164 = t149 * t154;
t153 = cos(pkin(7));
t163 = t150 * t153;
t151 = cos(pkin(12));
t162 = t151 * t153;
t152 = cos(pkin(11));
t161 = t152 * t154;
t147 = sin(pkin(12));
t143 = -t148 * t147 + t151 * t161;
t160 = -t143 * t153 + t152 * t165;
t145 = -t152 * t147 - t151 * t166;
t159 = t145 * t153 + t148 * t165;
t158 = cos(qJ(3));
t157 = cos(qJ(4));
t156 = sin(qJ(3));
t155 = sin(qJ(4));
t146 = -t147 * t166 + t152 * t151;
t144 = t147 * t161 + t148 * t151;
t142 = -t151 * t165 + t154 * t153;
t141 = -t145 * t149 + t148 * t163;
t140 = -t143 * t149 - t152 * t163;
t1 = [0, 0, t141, t146 * t156 - t159 * t158, 0 (t146 * t158 + t159 * t156) * t155 - t141 * t157; 0, 0, t140, t144 * t156 + t160 * t158, 0 (t144 * t158 - t160 * t156) * t155 - t140 * t157; 0, 0, t142, -t158 * t164 + (t147 * t156 - t158 * t162) * t150, 0 (t156 * t164 + (t147 * t158 + t156 * t162) * t150) * t155 - t142 * t157;];
Jg_rot  = t1;
