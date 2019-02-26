% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRP1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:41:34
% EndTime: 2019-02-26 19:41:34
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->20), mult. (98->46), div. (0->0), fcn. (140->12), ass. (0->28)
t149 = sin(pkin(11));
t155 = cos(pkin(6));
t167 = t149 * t155;
t150 = sin(pkin(7));
t151 = sin(pkin(6));
t166 = t150 * t151;
t165 = t150 * t155;
t154 = cos(pkin(7));
t164 = t151 * t154;
t152 = cos(pkin(12));
t163 = t152 * t154;
t153 = cos(pkin(11));
t162 = t153 * t155;
t148 = sin(pkin(12));
t144 = -t149 * t148 + t152 * t162;
t161 = -t144 * t154 + t153 * t166;
t146 = -t153 * t148 - t152 * t167;
t160 = t146 * t154 + t149 * t166;
t159 = cos(qJ(3));
t158 = cos(qJ(4));
t157 = sin(qJ(3));
t156 = sin(qJ(4));
t147 = -t148 * t167 + t153 * t152;
t145 = t148 * t162 + t149 * t152;
t143 = -t152 * t166 + t155 * t154;
t142 = -t146 * t150 + t149 * t164;
t141 = -t144 * t150 - t153 * t164;
t1 = [0, 0, t142, t147 * t157 - t160 * t159 (t147 * t159 + t160 * t157) * t156 - t142 * t158, 0; 0, 0, t141, t145 * t157 + t161 * t159 (t145 * t159 - t161 * t157) * t156 - t141 * t158, 0; 0, 0, t143, -t159 * t165 + (t148 * t157 - t159 * t163) * t151 (t157 * t165 + (t148 * t159 + t157 * t163) * t151) * t156 - t143 * t158, 0;];
Jg_rot  = t1;
